#!/usr/bin/env python3
"""Genera las figuras farmacogenómicas definitivas del TFM.

El script debe guardarse dentro de la carpeta ``scripts`` del repositorio.

Entradas:
    data_processed/master_her2_auc.csv
    results/tables/tables_compuesto/spearman_por_compuesto_completo.csv
    results/tables/tables_compuesto/mannwhitney_por_compuesto_completo.csv

Salidas:
    results/figures/figures_corrected/

Figuras generadas:
    Figura 2: distribución de ERBB2 y percentiles 20/80.
    Figura 3: relación entre ERBB2 y AUC de Dinaciclib.
    Figura 4: AUC de Dinaciclib en los grupos ERBB2-high y ERBB2-low.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


# -----------------------------------------------------------------------------
# 0. Rutas del proyecto
# -----------------------------------------------------------------------------
DIRECTORIO_PROYECTO = Path(__file__).resolve().parents[1]

ARCHIVO_MAESTRO = (
    DIRECTORIO_PROYECTO / "data_processed" / "master_her2_auc.csv"
)
DIRECTORIO_TABLAS = (
    DIRECTORIO_PROYECTO / "results" / "tables" / "tables_compuesto"
)
ARCHIVO_SPEARMAN = DIRECTORIO_TABLAS / "spearman_por_compuesto_completo.csv"
ARCHIVO_MANNWHITNEY = (
    DIRECTORIO_TABLAS / "mannwhitney_por_compuesto_completo.csv"
)
DIRECTORIO_FIGURAS = (
    DIRECTORIO_PROYECTO / "results" / "figures" / "figures_corrected"
)

VALORES_ESPERADOS = {
    "lineas": 22,
    "compuestos": 1383,
    "observaciones": 17253,
    "p20": 4.592383670461582,
    "p80": 7.386692396713945,
    "grupos": {"HER2_low": 5, "HER2_mid": 12, "HER2_high": 5},
    "dinaciclib_n": 21,
    "dinaciclib_rho": 0.778,
    "dinaciclib_fdr_spearman": 0.029,
}


# -----------------------------------------------------------------------------
# 1. Funciones auxiliares
# -----------------------------------------------------------------------------
def comprobar_archivo(ruta: Path) -> None:
    """Detiene la ejecución cuando falta un archivo requerido."""
    if not ruta.is_file():
        raise FileNotFoundError(f"No se encontró el archivo requerido:\n{ruta}")


def nombre_compuesto_corto(valor: object) -> str:
    """Elimina el identificador BRD y normaliza el nombre del compuesto."""
    return str(valor).split(" (BRD:", maxsplit=1)[0].strip().upper()


def guardar_figura(figura: plt.Figure, nombre: str) -> None:
    """Guarda cada figura en PNG a 300 dpi y en PDF vectorial."""
    DIRECTORIO_FIGURAS.mkdir(parents=True, exist_ok=True)
    figura.savefig(
        DIRECTORIO_FIGURAS / f"{nombre}.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    figura.savefig(
        DIRECTORIO_FIGURAS / f"{nombre}.pdf",
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(figura)


def añadir_marco(eje: plt.Axes) -> None:
    """Añade un marco fino alrededor del área de representación."""
    for borde in eje.spines.values():
        borde.set_visible(True)
        borde.set_color("#303030")
        borde.set_linewidth(0.9)


def obtener_fila_compuesto(
    tabla: pd.DataFrame,
    compuesto: str,
) -> pd.Series:
    """Localiza una única fila de un compuesto en una tabla de resultados."""
    columna_nombre = next(
        (
            columna
            for columna in ["Compuesto", "Compuesto completo", "compound"]
            if columna in tabla.columns
        ),
        None,
    )

    if columna_nombre is None:
        raise ValueError(
            "No se encontró una columna de compuesto en la tabla de resultados."
        )

    mascara = tabla[columna_nombre].map(nombre_compuesto_corto).eq(
        compuesto.upper()
    )

    if mascara.sum() != 1:
        raise ValueError(
            f"Se esperaba una fila para {compuesto}, pero se encontraron "
            f"{mascara.sum()}."
        )

    return tabla.loc[mascara].iloc[0]


def validar_dataset_maestro(
    datos: pd.DataFrame,
) -> tuple[pd.DataFrame, float, float]:
    """Comprueba que se está utilizando el dataset maestro definitivo."""
    columnas_necesarias = {
        "depmap_id",
        "compound",
        "AUC",
        "ERBB2_expr",
        "HER2_group",
    }
    columnas_faltantes = sorted(columnas_necesarias.difference(datos.columns))
    if columnas_faltantes:
        raise ValueError(
            "Faltan columnas en el dataset maestro: "
            + ", ".join(columnas_faltantes)
        )

    duplicados = datos.duplicated(["depmap_id", "compound"]).sum()
    if duplicados:
        raise ValueError(
            f"El dataset maestro contiene {duplicados} pares "
            "línea–compuesto duplicados."
        )

    asignaciones = datos[
        ["depmap_id", "ERBB2_expr", "HER2_group"]
    ].drop_duplicates()
    conflictos = asignaciones.groupby("depmap_id").size()
    conflictos = conflictos[conflictos != 1]
    if not conflictos.empty:
        raise ValueError(
            "Existen líneas con más de un valor de ERBB2 o grupo asignado: "
            + ", ".join(conflictos.index.astype(str))
        )

    lineas = asignaciones.sort_values("ERBB2_expr").reset_index(drop=True)
    p20 = float(lineas["ERBB2_expr"].quantile(0.20))
    p80 = float(lineas["ERBB2_expr"].quantile(0.80))

    observados = {
        "lineas": lineas["depmap_id"].nunique(),
        "compuestos": datos["compound"].nunique(),
        "observaciones": len(datos),
    }
    problemas: list[str] = []

    for metrica, observado in observados.items():
        esperado = VALORES_ESPERADOS[metrica]
        if observado != esperado:
            problemas.append(
                f"{metrica}={observado}; valor esperado={esperado}"
            )

    if not np.isclose(p20, VALORES_ESPERADOS["p20"], atol=1e-10):
        problemas.append(
            f"P20={p20:.12f}; valor esperado={VALORES_ESPERADOS['p20']:.12f}"
        )
    if not np.isclose(p80, VALORES_ESPERADOS["p80"], atol=1e-10):
        problemas.append(
            f"P80={p80:.12f}; valor esperado={VALORES_ESPERADOS['p80']:.12f}"
        )

    grupos = lineas["HER2_group"].value_counts().to_dict()
    if grupos != VALORES_ESPERADOS["grupos"]:
        problemas.append(
            f"grupos={grupos}; valor esperado={VALORES_ESPERADOS['grupos']}"
        )

    if problemas:
        raise RuntimeError(
            "El archivo no coincide con el dataset maestro definitivo:\n- "
            + "\n- ".join(problemas)
        )

    return lineas, p20, p80


def validar_resultados_dinaciclib(
    datos_dinaciclib: pd.DataFrame,
    resultados_spearman: pd.DataFrame,
    resultados_mw: pd.DataFrame,
) -> tuple[pd.Series, pd.Series]:
    """Comprueba la concordancia de las tablas con los datos de Dinaciclib."""
    if len(datos_dinaciclib) != VALORES_ESPERADOS["dinaciclib_n"]:
        raise RuntimeError(
            f"Dinaciclib contiene {len(datos_dinaciclib)} observaciones; "
            f"se esperaban {VALORES_ESPERADOS['dinaciclib_n']}."
        )

    if datos_dinaciclib["depmap_id"].nunique() != len(datos_dinaciclib):
        raise RuntimeError("Dinaciclib contiene líneas celulares duplicadas.")

    fila_spearman = obtener_fila_compuesto(resultados_spearman, "DINACICLIB")
    fila_mw = obtener_fila_compuesto(resultados_mw, "DINACICLIB")

    columnas_spearman = {"Rho de Spearman", "p-valor", "FDR"}
    faltantes_spearman = columnas_spearman.difference(resultados_spearman.columns)
    if faltantes_spearman:
        raise ValueError(
            "Faltan columnas en la tabla de Spearman: "
            + ", ".join(sorted(faltantes_spearman))
        )

    columnas_mw = {"p-valor", "FDR"}
    faltantes_mw = columnas_mw.difference(resultados_mw.columns)
    if faltantes_mw:
        raise ValueError(
            "Faltan columnas en la tabla de Mann–Whitney: "
            + ", ".join(sorted(faltantes_mw))
        )

    rho_calculado, p_spearman_calculado = stats.spearmanr(
        datos_dinaciclib["ERBB2_expr"],
        datos_dinaciclib["AUC"],
    )
    if not np.isclose(float(fila_spearman["Rho de Spearman"]), rho_calculado):
        raise RuntimeError(
            "El rho de Spearman de la tabla no coincide con los datos de Dinaciclib."
        )
    if not np.isclose(float(fila_spearman["p-valor"]), p_spearman_calculado):
        raise RuntimeError(
            "El p-valor de Spearman no coincide con los datos de Dinaciclib."
        )
    if not np.isclose(
        float(fila_spearman["Rho de Spearman"]),
        VALORES_ESPERADOS["dinaciclib_rho"],
        atol=0.001,
    ):
        raise RuntimeError("El rho de Dinaciclib no coincide con el resultado final.")
    if not np.isclose(
        float(fila_spearman["FDR"]),
        VALORES_ESPERADOS["dinaciclib_fdr_spearman"],
        atol=0.001,
    ):
        raise RuntimeError("El FDR de Dinaciclib no coincide con el resultado final.")

    extremos = datos_dinaciclib.loc[
        datos_dinaciclib["HER2_group"].isin(["HER2_high", "HER2_low"])
    ]
    grupo_high = extremos.loc[extremos["HER2_group"].eq("HER2_high"), "AUC"]
    grupo_low = extremos.loc[extremos["HER2_group"].eq("HER2_low"), "AUC"]
    p_mw_calculado = stats.mannwhitneyu(
        grupo_high,
        grupo_low,
        alternative="two-sided",
    ).pvalue
    if not np.isclose(float(fila_mw["p-valor"]), p_mw_calculado):
        raise RuntimeError(
            "El p-valor de Mann–Whitney no coincide con los datos de Dinaciclib."
        )

    return fila_spearman, fila_mw


# -----------------------------------------------------------------------------
# 2. Figuras
# -----------------------------------------------------------------------------
def crear_distribucion_erbb2(
    lineas: pd.DataFrame,
    p20: float,
    p80: float,
) -> None:
    """Genera la figura 2: distribución de ERBB2 y grupos extremos."""
    datos_figura = lineas.copy()
    etiquetas = {
        "HER2_low": "≤ P20: expresión baja (n = 5)",
        "HER2_mid": "Entre P20 y P80 (n = 12)",
        "HER2_high": "≥ P80: expresión alta (n = 5)",
    }
    orden = [
        "≤ P20: expresión baja (n = 5)",
        "Entre P20 y P80 (n = 12)",
        "≥ P80: expresión alta (n = 5)",
    ]
    paleta = {
        orden[0]: "#4C78A8",
        orden[1]: "#9AA0A6",
        orden[2]: "#D65F5F",
    }
    datos_figura["Grupo"] = datos_figura["HER2_group"].map(etiquetas)

    figura, eje = plt.subplots(figsize=(8.4, 4.6))
    sns.swarmplot(
        data=datos_figura,
        x="ERBB2_expr",
        y="Grupo",
        hue="Grupo",
        order=orden,
        hue_order=orden,
        palette=paleta,
        size=8,
        edgecolor="white",
        linewidth=0.7,
        legend=False,
        ax=eje,
    )
    eje.axvline(
        p20,
        color="#2647A0",
        linestyle="--",
        linewidth=1.8,
        label=f"Percentil 20 = {p20:.3f}",
    )
    eje.axvline(
        p80,
        color="#2A8F3A",
        linestyle="--",
        linewidth=1.8,
        label=f"Percentil 80 = {p80:.3f}",
    )
    eje.set_title(
        "Distribución de la expresión de ERBB2 en las 22 líneas celulares",
        pad=14,
    )
    eje.set_xlabel("Expresión de ERBB2, log2(TPM+1)")
    eje.set_ylabel("Clasificación según percentiles")
    eje.legend(frameon=False)
    eje.grid(axis="x", alpha=0.2)
    eje.grid(axis="y", visible=False)
    añadir_marco(eje)
    figura.tight_layout()
    guardar_figura(
        figura,
        "figura_02_distribucion_ERBB2_percentiles_FINAL",
    )


def crear_scatter_dinaciclib(
    datos_dinaciclib: pd.DataFrame,
    fila_spearman: pd.Series,
) -> None:
    """Genera la figura 3: relación entre ERBB2 y AUC de Dinaciclib."""
    rho = float(fila_spearman["Rho de Spearman"])
    p_valor = float(fila_spearman["p-valor"])
    fdr = float(fila_spearman["FDR"])

    figura, eje = plt.subplots(figsize=(7.4, 5.1))
    sns.regplot(
        data=datos_dinaciclib,
        x="ERBB2_expr",
        y="AUC",
        scatter_kws={
            "s": 48,
            "alpha": 0.85,
            "color": "#4C88B7",
            "edgecolor": "white",
        },
        line_kws={"color": "#173B73", "linewidth": 2},
        ci=95,
        seed=20260903,
        ax=eje,
    )
    eje.set_title("Relación entre la expresión de ERBB2 y el AUC de Dinaciclib")
    eje.set_xlabel("Expresión de ERBB2, log2(TPM+1)")
    eje.set_ylabel("AUC de PRISM")
    eje.text(
        0.04,
        0.96,
        f"n = {len(datos_dinaciclib)}\n"
        f"ρ de Spearman = {rho:.3f}\n"
        f"p = {p_valor:.2e}\n"
        f"FDR = {fdr:.3f}",
        transform=eje.transAxes,
        va="top",
        ha="left",
        fontsize=10.5,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "alpha": 0.85,
            "edgecolor": "none",
        },
    )
    eje.text(
        0.98,
        0.03,
        "Ajuste lineal descriptivo",
        transform=eje.transAxes,
        ha="right",
        va="bottom",
        fontsize=8.5,
        color="#555555",
    )
    eje.grid(alpha=0.18)
    añadir_marco(eje)
    figura.tight_layout()
    guardar_figura(figura, "figura_03_scatter_ERBB2_Dinaciclib")


def crear_boxplot_dinaciclib(
    datos_dinaciclib: pd.DataFrame,
    fila_mw: pd.Series,
) -> None:
    """Genera la figura 4: AUC de Dinaciclib en los grupos extremos."""
    extremos = datos_dinaciclib.loc[
        datos_dinaciclib["HER2_group"].isin(["HER2_high", "HER2_low"])
    ].copy()
    extremos["Grupo"] = extremos["HER2_group"].map(
        {"HER2_high": "ERBB2-high", "HER2_low": "ERBB2-low"}
    )

    grupo_high = extremos.loc[extremos["Grupo"].eq("ERBB2-high"), "AUC"]
    grupo_low = extremos.loc[extremos["Grupo"].eq("ERBB2-low"), "AUC"]
    orden = [
        f"ERBB2-high\n(n = {len(grupo_high)})",
        f"ERBB2-low\n(n = {len(grupo_low)})",
    ]
    extremos["Grupo_n"] = extremos["Grupo"].map(
        {"ERBB2-high": orden[0], "ERBB2-low": orden[1]}
    )

    np.random.seed(20260903)
    figura, eje = plt.subplots(figsize=(7.2, 5.0))
    sns.boxplot(
        data=extremos,
        x="Grupo_n",
        y="AUC",
        order=orden,
        color="#A9CEE0",
        width=0.55,
        showfliers=False,
        linewidth=1.2,
        ax=eje,
    )
    sns.stripplot(
        data=extremos,
        x="Grupo_n",
        y="AUC",
        order=orden,
        color="#3E7FA6",
        edgecolor="white",
        linewidth=0.6,
        alpha=0.9,
        jitter=0.08,
        size=7,
        ax=eje,
    )
    figura.suptitle(
        "AUC de Dinaciclib por grupo de expresión de ERBB2",
        fontsize=14,
        y=0.975,
    )
    eje.set_xlabel("Grupo transcriptómico")
    eje.set_ylabel("AUC de PRISM")
    eje.grid(axis="y", alpha=0.18)
    añadir_marco(eje)
    figura.tight_layout(rect=[0, 0, 1, 0.91])
    guardar_figura(
        figura,
        "figura_04_boxplot_Dinaciclib_ERBB2_FINAL",
    )

    print(
        "Dinaciclib (Mann–Whitney bilateral): "
        f"p = {float(fila_mw['p-valor']):.4f}; "
        f"FDR = {float(fila_mw['FDR']):.3f}"
    )


# -----------------------------------------------------------------------------
# 3. Flujo principal
# -----------------------------------------------------------------------------
def main() -> None:
    """Carga y valida las entradas y genera las figuras definitivas."""
    comprobar_archivo(ARCHIVO_MAESTRO)
    comprobar_archivo(ARCHIVO_SPEARMAN)
    comprobar_archivo(ARCHIVO_MANNWHITNEY)

    datos = pd.read_csv(ARCHIVO_MAESTRO, low_memory=False)
    datos["ERBB2_expr"] = pd.to_numeric(datos["ERBB2_expr"], errors="coerce")
    datos["AUC"] = pd.to_numeric(datos["AUC"], errors="coerce")
    datos = datos.dropna(
        subset=["depmap_id", "compound", "ERBB2_expr", "AUC", "HER2_group"]
    )

    lineas, p20, p80 = validar_dataset_maestro(datos)

    datos["compound_short"] = datos["compound"].map(nombre_compuesto_corto)
    datos_dinaciclib = datos.loc[
        datos["compound_short"].eq("DINACICLIB")
    ].copy()

    resultados_spearman = pd.read_csv(ARCHIVO_SPEARMAN)
    resultados_mw = pd.read_csv(ARCHIVO_MANNWHITNEY)
    fila_spearman, fila_mw = validar_resultados_dinaciclib(
        datos_dinaciclib,
        resultados_spearman,
        resultados_mw,
    )

    sns.set_theme(style="whitegrid", context="notebook")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "xtick.labelsize": 10.5,
            "ytick.labelsize": 10.5,
        }
    )

    crear_distribucion_erbb2(lineas, p20, p80)
    crear_scatter_dinaciclib(datos_dinaciclib, fila_spearman)
    crear_boxplot_dinaciclib(datos_dinaciclib, fila_mw)

    print("\nValidación correcta:")
    print(f"- Líneas celulares: {lineas['depmap_id'].nunique()}")
    print(f"- Compuestos: {datos['compound'].nunique()}")
    print(f"- Observaciones: {len(datos)}")
    print(f"- P20: {p20:.12f}")
    print(f"- P80: {p80:.12f}")
    print(f"- Figuras guardadas en: {DIRECTORIO_FIGURAS}")


if __name__ == "__main__":
    main()
