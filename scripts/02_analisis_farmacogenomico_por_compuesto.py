#!/usr/bin/env python3
"""Este script realiza el análisis farmacogenómico por compuesto y calcula las correlaciones
de Spearman, las comparaciones de Mann–Whitney, el ajuste de Benjamini–Hochberg
y las tablas empleadas en el documento."""

# Cargar librerias y paquetes
from math import isclose
from pathlib import Path
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests


# 0. Rutas y parámetros
DIRECTORIO_PROYECTO = Path(os.environ.get("TFM_DIR", Path(__file__).resolve().parents[1])).expanduser().resolve()
# Todas las salidas van a una copia de comprobación; el repositorio queda intacto.
DIRECTORIO_SALIDAS = DIRECTORIO_PROYECTO / "comprobacion_reproducibilidad"
ARCHIVO_MAESTRO = (
    DIRECTORIO_SALIDAS / "datos_procesados" / "master_her2_auc.csv"
)
DIRECTORIO_TABLAS = (
    DIRECTORIO_SALIDAS / "resultados" / "farmacogenomica"
)
DIRECTORIO_APOYO = DIRECTORIO_SALIDAS / "resultados" / "tablas_apoyo"
DIRECTORIO_CONTROL = DIRECTORIO_SALIDAS / "resultados" / "control_calidad"
DIRECTORIO_TABLAS_TFM = DIRECTORIO_SALIDAS / "tablas"

DIRECTORIO_FIGURAS = (
    DIRECTORIO_SALIDAS / "figuras" / "complementarias"
)

MINIMO_SPEARMAN = 10
MINIMO_POR_GRUPO = 3

VALORES_ESPERADOS = {
    "lineas": 22,
    "compuestos": 1383,
    "observaciones": 17253,
    "spearman_evaluados": 883,
    "mannwhitney_evaluados": 572,
    "spearman_fdr_005": 1,
    "mannwhitney_fdr_005": 0,
}

NOMBRES_COMPUESTOS_ES = {
    "DINACICLIB": "Dinaciclib",
    "ALPELISIB": "Alpelisib",
    "VELBAN": "Velban",
    "EMETINE": "Emetina",
    "CANERTINIB": "Canertinib",
    "CLONAZEPAM": "Clonazepam",
    "POZIOTINIB": "Poziotinib",
    "BETULINIC-ACID": "Ácido betulínico",
    "LAPATINIB": "Lapatinib",
    "TUCATINIB": "Tucatinib",
    "GEFITINIB": "Gefitinib",
}

MECANISMOS_ACCION = {
    "DINACICLIB": "Inhibidor de CDK",
    "ALPELISIB": "Inhibidor de PI3Kα",
    "VELBAN": "Inhibidor de microtúbulos",
    "EMETINE": "Inhibidor de la síntesis proteica",
    "CANERTINIB": "Inhibidor pan-HER/EGFR",
    "CLONAZEPAM": "Modulador del receptor GABA-A",
    "POZIOTINIB": "Inhibidor de EGFR/HER2",
    "BETULINIC-ACID": "Triterpenoide con actividad antitumoral",
    "LAPATINIB": "Inhibidor de HER2/EGFR",
    "TUCATINIB": "Inhibidor selectivo de HER2",
    "GEFITINIB": "Inhibidor de EGFR",
}

COMPUESTOS_INTERES_MW = [
    "DINACICLIB",
    "ALPELISIB",
    "CANERTINIB",
    "GEFITINIB",
    "LAPATINIB",
    "TUCATINIB",
    "POZIOTINIB",
]


# 1. Funciones auxiliares
def nombre_compuesto_corto(nombre: object) -> str:
    """Elimina el identificador BRD y normaliza el nombre del compuesto."""
    return str(nombre).split(" (BRD:", maxsplit=1)[0].strip().upper()


def etiqueta_grupo_her2(valor: str) -> str:
    """Convierte las etiquetas internas en etiquetas para las figuras."""
    etiquetas = {
        "HER2_high": "HER2-high",
        "HER2_low": "HER2-low",
        "HER2_mid": "HER2-intermedio",
    }
    return etiquetas.get(valor, valor)


def guardar_tabla_tfm(tabla: pd.DataFrame, ruta_base: Path) -> None:
    """Guarda CSV con punto y coma y coma decimal para Excel/Numbers."""
    ruta_base.parent.mkdir(parents=True, exist_ok=True)
    tabla.to_csv(ruta_base.with_suffix(".csv"), index=False,
                 sep=";", decimal=",", encoding="utf-8-sig")


def validar_dataset_maestro(datos: pd.DataFrame) -> None:
    """Comprueba que se está utilizando el dataset maestro definitivo."""
    columnas_necesarias = {
        "depmap_id",
        "compound",
        "AUC",
        "ERBB2_expr",
        "HER2_group",
    }
    faltantes = sorted(columnas_necesarias.difference(datos.columns))
    if faltantes:
        raise ValueError(
            "Faltan columnas en el dataset maestro: " + ", ".join(faltantes)
        )

    duplicados = datos.duplicated(["depmap_id", "compound"]).sum()
    observados = {
        "lineas": datos["depmap_id"].nunique(),
        "compuestos": datos["compound"].nunique(),
        "observaciones": len(datos),
    }

    errores = []
    for metrica, valor in observados.items():
        if valor != VALORES_ESPERADOS[metrica]:
            errores.append(
                f"{metrica}: {valor} (esperado: {VALORES_ESPERADOS[metrica]})"
            )

    if duplicados:
        errores.append(f"pares línea–compuesto duplicados: {duplicados}")

    if not np.isfinite(datos[["AUC", "ERBB2_expr"]].to_numpy()).all():
        errores.append("valores no finitos en AUC o ERBB2")
    asignaciones = datos[["depmap_id", "ERBB2_expr", "HER2_group"]].drop_duplicates()
    if asignaciones["depmap_id"].duplicated().any():
        errores.append("una línea tiene múltiples expresiones o grupos")
    lineas = datos[
        ["depmap_id", "ERBB2_expr", "HER2_group"]
    ].drop_duplicates(subset="depmap_id")
    grupos = lineas["HER2_group"].value_counts().to_dict()
    grupos_esperados = {"HER2_low": 5, "HER2_mid": 12, "HER2_high": 5}
    if grupos != grupos_esperados:
        errores.append(
            f"distribución de grupos: {grupos} (esperada: {grupos_esperados})"
        )

    p20 = lineas["ERBB2_expr"].quantile(0.20)
    p80 = lineas["ERBB2_expr"].quantile(0.80)
    if not isclose(p20, 4.592383670461582, abs_tol=1e-10):
        errores.append(f"P20 inesperado: {p20}")
    if not isclose(p80, 7.386692396713945, abs_tol=1e-10):
        errores.append(f"P80 inesperado: {p80}")

    grupos_calculados = np.where(lineas["ERBB2_expr"] <= p20, "HER2_low",
                                 np.where(lineas["ERBB2_expr"] >= p80, "HER2_high", "HER2_mid"))
    if not np.array_equal(grupos_calculados, lineas["HER2_group"].to_numpy()):
        errores.append("grupos incompatibles con los percentiles")
    if errores:
        raise RuntimeError(
            "El archivo no coincide con el dataset maestro definitivo:\n- "
            + "\n- ".join(errores)
        )


def calcular_spearman(datos: pd.DataFrame) -> pd.DataFrame:
    """Calcula Spearman por compuesto y ajusta los p-valores mediante BH."""
    resultados = []

    for compuesto, subdatos in datos.groupby("compound"):
        subdatos = subdatos.dropna(subset=["ERBB2_expr", "AUC"])
        n_total = len(subdatos)
        n_lineas = subdatos["depmap_id"].nunique()

        if n_lineas < MINIMO_SPEARMAN:
            continue

        rho, p_valor = stats.spearmanr(
            subdatos["ERBB2_expr"], subdatos["AUC"]
        )
        if pd.isna(rho) or pd.isna(p_valor):
            continue

        resultados.append(
            {
                "Compuesto completo": compuesto,
                "Compuesto": nombre_compuesto_corto(compuesto),
                "n": int(n_total),
                "Líneas celulares": int(n_lineas),
                "Rho de Spearman": rho,
                "p-valor": p_valor,
            }
        )

    tabla = pd.DataFrame(resultados)
    if tabla.empty:
        raise RuntimeError("No se obtuvieron resultados válidos de Spearman.")

    tabla["FDR"] = multipletests(tabla["p-valor"], method="fdr_bh")[1]
    return tabla.sort_values(
        ["FDR", "p-valor", "Compuesto"]
    ).reset_index(drop=True)


def calcular_mannwhitney(
    datos: pd.DataFrame,
    aplicar_minimo: bool = True,
) -> pd.DataFrame:
    """Compara HER2-high y HER2-low individualmente para cada compuesto."""
    resultados = []

    for compuesto, subdatos in datos.groupby("compound"):
        grupo_high = subdatos.loc[
            subdatos["HER2_group"].eq("HER2_high"), "AUC"
        ].dropna()
        grupo_low = subdatos.loc[
            subdatos["HER2_group"].eq("HER2_low"), "AUC"
        ].dropna()

        n_high = len(grupo_high)
        n_low = len(grupo_low)

        if aplicar_minimo:
            if n_high < MINIMO_POR_GRUPO or n_low < MINIMO_POR_GRUPO:
                continue
        elif n_high == 0 or n_low == 0:
            continue

        _, p_valor = stats.mannwhitneyu(
            grupo_high,
            grupo_low,
            alternative="two-sided",
            method="auto",
        )

        resultados.append(
            {
                "Compuesto completo": compuesto,
                "Compuesto": nombre_compuesto_corto(compuesto),
                "n HER2-high": int(n_high),
                "n HER2-low": int(n_low),
                "Mediana AUC HER2-high": grupo_high.median(),
                "Mediana AUC HER2-low": grupo_low.median(),
                "Diferencia de medianas": (
                    grupo_high.median() - grupo_low.median()
                ),
                "p-valor": p_valor,
                "Criterio n>=3": (
                    "Sí"
                    if n_high >= MINIMO_POR_GRUPO
                    and n_low >= MINIMO_POR_GRUPO
                    else "No"
                ),
            }
        )

    tabla = pd.DataFrame(resultados)
    if tabla.empty:
        raise RuntimeError("No se obtuvieron resultados válidos de Mann–Whitney.")

    tabla["FDR"] = multipletests(tabla["p-valor"], method="fdr_bh")[1]
    return tabla.sort_values(
        ["p-valor", "FDR", "Compuesto"]
    ).reset_index(drop=True)


def validar_resultados(
    spearman: pd.DataFrame,
    mannwhitney: pd.DataFrame,
) -> None:
    """Valida los resultados fundamentales informados en el TFM."""
    errores = []
    comprobaciones = {
        "spearman_evaluados": len(spearman),
        "mannwhitney_evaluados": len(mannwhitney),
        "spearman_fdr_005": int((spearman["FDR"] < 0.05).sum()),
        "mannwhitney_fdr_005": int((mannwhitney["FDR"] < 0.05).sum()),
    }

    for metrica, observado in comprobaciones.items():
        esperado = VALORES_ESPERADOS[metrica]
        if observado != esperado:
            errores.append(f"{metrica}: {observado} (esperado: {esperado})")

    dinaciclib = spearman.loc[spearman["Compuesto"].eq("DINACICLIB")]
    if len(dinaciclib) != 1:
        errores.append("no se encontró una única fila para Dinaciclib")
    else:
        fila = dinaciclib.iloc[0]
        if not isclose(float(fila["Rho de Spearman"]), 0.778, abs_tol=0.001):
            errores.append(
                f"rho de Dinaciclib inesperado: {fila['Rho de Spearman']}"
            )
        if not isclose(float(fila["FDR"]), 0.029, abs_tol=0.001):
            errores.append(f"FDR de Dinaciclib inesperado: {fila['FDR']}")

    if errores:
        raise RuntimeError(
            "Los resultados no coinciden con la ejecución definitiva:\n- "
            + "\n- ".join(errores)
        )


def crear_figura_compuestos(datos_extremos: pd.DataFrame) -> None:
    """Genera la comparación final de AUC para cinco compuestos de interés."""
    orden_original = [
        "CANERTINIB",
        "GEFITINIB",
        "TUCATINIB",
        "POZIOTINIB",
        "DINACICLIB",
    ]
    orden_es = [NOMBRES_COMPUESTOS_ES[x] for x in orden_original]

    datos_figura = datos_extremos.loc[
        datos_extremos["compound"]
        .map(nombre_compuesto_corto)
        .isin(orden_original)
    ].copy()
    datos_figura["Compuesto"] = (
        datos_figura["compound"]
        .map(nombre_compuesto_corto)
        .map(NOMBRES_COMPUESTOS_ES)
    )
    datos_figura["Grupo HER2"] = datos_figura["HER2_group"].map(
        etiqueta_grupo_her2
    )

    encontrados = set(datos_figura["Compuesto"].dropna())
    if encontrados != set(orden_es):
        raise RuntimeError(
            "No se encontraron todos los compuestos necesarios para la figura 5."
        )

    paleta = {"HER2-high": "#2B6CB0", "HER2-low": "#90CDF4"}
    np.random.seed(20260903)
    figura, eje = plt.subplots(figsize=(8, 5.5))

    sns.boxplot(
        data=datos_figura,
        y="Compuesto",
        x="AUC",
        hue="Grupo HER2",
        order=orden_es,
        hue_order=["HER2-high", "HER2-low"],
        width=0.65,
        fliersize=0,
        palette=paleta,
        ax=eje,
    )
    sns.stripplot(
        data=datos_figura,
        y="Compuesto",
        x="AUC",
        hue="Grupo HER2",
        order=orden_es,
        hue_order=["HER2-high", "HER2-low"],
        dodge=True,
        jitter=0.18,
        alpha=0.75,
        size=4,
        edgecolor="black",
        linewidth=0.4,
        palette=paleta,
        ax=eje,
    )

    eje.xaxis.grid(
        True,
        linestyle=":",
        color="lightgray",
        linewidth=0.8,
        alpha=0.7,
    )
    for posicion in range(len(orden_es) - 1):
        eje.axhline(
            posicion + 0.5,
            color="lightgray",
            linestyle="--",
            linewidth=0.7,
            alpha=0.5,
            zorder=0,
        )

    eje.set_title(
        "Comparación de la sensibilidad farmacológica entre grupos HER2"
    )
    eje.set_xlabel("Sensibilidad farmacológica (AUC)")
    eje.set_ylabel("Compuesto")
    eje.tick_params(axis="x", rotation=45)

    manejadores, etiquetas = eje.get_legend_handles_labels()
    eje.legend(
        manejadores[:2],
        etiquetas[:2],
        title="Grupo HER2",
        frameon=False,
    )

    figura.tight_layout()
    ruta_base = DIRECTORIO_FIGURAS / "comparacion_compuestos_interes_HER2"
    figura.savefig(ruta_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    figura.savefig(ruta_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figura)


# 2. Flujo principal
def main() -> None:
    """Ejecuta el análisis completo por compuesto."""
    if not ARCHIVO_MAESTRO.is_file():
        raise FileNotFoundError(
            f"No se encontró el dataset maestro requerido:\n{ARCHIVO_MAESTRO}"
        )

    for carpeta in (DIRECTORIO_TABLAS, DIRECTORIO_APOYO, DIRECTORIO_CONTROL, DIRECTORIO_TABLAS_TFM):
        carpeta.mkdir(parents=True, exist_ok=True)
    DIRECTORIO_FIGURAS.mkdir(parents=True, exist_ok=True)

    datos = pd.read_csv(ARCHIVO_MAESTRO, low_memory=False)
    validar_dataset_maestro(datos)
    datos["ERBB2_expr"] = pd.to_numeric(datos["ERBB2_expr"], errors="coerce")
    datos["AUC"] = pd.to_numeric(datos["AUC"], errors="coerce")
    datos = datos.dropna(
        subset=["depmap_id", "compound", "ERBB2_expr", "AUC", "HER2_group"]
    )
    datos = datos.loc[
        datos["HER2_group"].isin(["HER2_high", "HER2_low", "HER2_mid"])
    ].copy()
    validar_dataset_maestro(datos)

    # Correlaciones de Spearman por compuesto.
    spearman_completo = calcular_spearman(datos)
    datos_extremos = datos.loc[datos["HER2_group"].isin(["HER2_high", "HER2_low"])].copy()
    mw_completo = calcular_mannwhitney(datos_extremos, aplicar_minimo=True)
    validar_resultados(spearman_completo, mw_completo)
    spearman_completo.to_csv(
        DIRECTORIO_TABLAS / "spearman_por_compuesto.csv",
        index=False,
    )

    tabla_2 = spearman_completo.head(10).copy()
    compuestos_originales = tabla_2["Compuesto"].copy()
    tabla_2["Compuesto"] = (
        compuestos_originales.map(NOMBRES_COMPUESTOS_ES)
        .fillna(compuestos_originales)
    )
    tabla_2["Mecanismo de acción"] = compuestos_originales.map(
        MECANISMOS_ACCION
    )
    tabla_2 = tabla_2[
        [
            "Compuesto",
            "Mecanismo de acción",
            "Líneas celulares",
            "Rho de Spearman",
            "p-valor",
            "FDR",
        ]
    ]
    tabla_2["Rho de Spearman"] = tabla_2["Rho de Spearman"].round(3)
    tabla_2["p-valor"] = tabla_2["p-valor"].round(6)
    tabla_2["FDR"] = tabla_2["FDR"].round(3)
    guardar_tabla_tfm(
        tabla_2,
        DIRECTORIO_TABLAS_TFM / "tabla_03_asociaciones_spearman",
    )

    # Mann–Whitney principal calculado y validado antes de exportar.
    mw_completo.to_csv(
        DIRECTORIO_TABLAS / "mannwhitney_por_compuesto.csv",
        index=False,
    )

    # Tabla exploratoria de los compuestos discutidos en el TFM. El FDR procede
    # siempre del análisis global, nunca de este subconjunto seleccionado.
    datos_interes = datos_extremos.loc[
        datos_extremos["compound"]
        .map(nombre_compuesto_corto)
        .isin(COMPUESTOS_INTERES_MW)
    ].copy()
    mw_interes = calcular_mannwhitney(datos_interes, aplicar_minimo=False)
    mw_interes = mw_interes.drop(columns="FDR")
    mw_interes = mw_interes.merge(
        mw_completo[["Compuesto completo", "FDR"]],
        on="Compuesto completo",
        how="left",
        validate="one_to_one",
    ).rename(columns={"FDR": "FDR (global)"})
    mw_interes.to_csv(
        DIRECTORIO_APOYO / "mannwhitney_compuestos_interes_exploratorio.csv",
        index=False,
    )

    tabla_3 = mw_interes.copy()
    compuestos_originales = tabla_3["Compuesto"].copy()
    tabla_3["Compuesto"] = (
        compuestos_originales.map(NOMBRES_COMPUESTOS_ES)
        .fillna(compuestos_originales)
    )
    tabla_3["Mecanismo de acción"] = compuestos_originales.map(
        MECANISMOS_ACCION
    )
    tabla_3 = tabla_3[
        [
            "Compuesto",
            "Mecanismo de acción",
            "n HER2-high",
            "n HER2-low",
            "Mediana AUC HER2-high",
            "Mediana AUC HER2-low",
            "Diferencia de medianas",
            "p-valor",
            "FDR (global)",
            "Criterio n>=3",
        ]
    ]
    columnas_redondear = [
        "Mediana AUC HER2-high",
        "Mediana AUC HER2-low",
        "Diferencia de medianas",
        "p-valor",
        "FDR (global)",
    ]
    tabla_3[columnas_redondear] = tabla_3[columnas_redondear].round(4)
    guardar_tabla_tfm(
        tabla_3,
        DIRECTORIO_TABLAS_TFM / "tabla_04_comparacion_mannwhitney",
    )

    diagnostico_mw = pd.DataFrame(
        {
            "Métrica": [
                "Compuestos analizados",
                "p-valores únicos",
                "FDR únicos",
                "p < 0,05",
                "FDR < 0,05",
                "FDR < 0,10",
            ],
            "Valor": [
                len(mw_completo),
                mw_completo["p-valor"].nunique(),
                mw_completo["FDR"].nunique(),
                int((mw_completo["p-valor"] < 0.05).sum()),
                int((mw_completo["FDR"] < 0.05).sum()),
                int((mw_completo["FDR"] < 0.10).sum()),
            ],
        }
    )
    diagnostico_mw.to_csv(
        DIRECTORIO_CONTROL / "diagnostico_pvalores_FDR_mannwhitney.csv",
        index=False,
    )

    # Integración de Spearman y Mann–Whitney para la interpretación biológica.
    tabla_integrada = spearman_completo.merge(
        mw_completo,
        on=["Compuesto completo", "Compuesto"],
        how="outer",
        suffixes=("_Spearman", "_MannWhitney"),
        validate="one_to_one",
    )
    tabla_integrada["Mecanismo de acción"] = tabla_integrada["Compuesto"].map(
        MECANISMOS_ACCION
    )
    tabla_integrada.to_csv(
        DIRECTORIO_TABLAS / "resultados_integrados_por_compuesto.csv",
        index=False,
    )

    tabla_biologica = tabla_integrada.loc[
        tabla_integrada["Compuesto"].isin(COMPUESTOS_INTERES_MW)
    ].copy()
    columnas_biologicas = [
        "Compuesto",
        "Mecanismo de acción",
        "n",
        "Rho de Spearman",
        "p-valor_Spearman",
        "FDR_Spearman",
        "n HER2-high",
        "n HER2-low",
        "Diferencia de medianas",
        "p-valor_MannWhitney",
        "FDR_MannWhitney",
    ]
    columnas_biologicas = [
        columna
        for columna in columnas_biologicas
        if columna in tabla_biologica.columns
    ]
    tabla_biologica = tabla_biologica[columnas_biologicas]
    for columna in tabla_biologica.select_dtypes(include="number").columns:
        tabla_biologica[columna] = tabla_biologica[columna].round(4)
    tabla_biologica.to_csv(
        DIRECTORIO_APOYO / "tabla_apoyo_interpretacion_biologica.csv",
        index=False,
    )

    sns.set_theme(style="whitegrid", context="notebook")
    crear_figura_compuestos(datos_extremos)

    print("\nAnálisis por compuesto completado y validado:")
    print(f"- Compuestos evaluados mediante Spearman: {len(spearman_completo)}")
    print(f"- Spearman con FDR < 0,05: {(spearman_completo['FDR'] < 0.05).sum()}")
    print(f"- Compuestos evaluados mediante Mann–Whitney: {len(mw_completo)}")
    print(f"- Mann–Whitney con FDR < 0,05: {(mw_completo['FDR'] < 0.05).sum()}")
    print(f"- Tablas guardadas en: {DIRECTORIO_TABLAS}")
    print(f"- Figura guardada en: {DIRECTORIO_FIGURAS}")


if __name__ == "__main__":
    main()
