#!/usr/bin/env python3
"""Construye el conjunto de datos maestro ERBB2–AUC del TFM."""

# Cargar paquetes necesarios
from math import isclose
from pathlib import Path
import os
import pandas as pd


# 0. Rutas del proyecto
DIRECTORIO_PROYECTO = Path(os.environ.get("TFM_DIR", Path(__file__).resolve().parents[1])).expanduser().resolve()
# Todas las salidas van a una copia de comprobación; el repositorio queda intacto.
DIRECTORIO_SALIDAS = DIRECTORIO_PROYECTO / "comprobacion_reproducibilidad"

ARCHIVO_AUC = (
    DIRECTORIO_PROYECTO
    / "datos_originales"
    / "Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_subsetted.csv"
)
ARCHIVO_EXPRESION = (
    DIRECTORIO_PROYECTO
    / "datos_originales"
    / "Expression_Public_25Q3_subsetted.csv"
)

ARCHIVO_MAESTRO = (
    DIRECTORIO_SALIDAS / "datos_procesados" / "master_her2_auc.csv"
)
TABLA_RESUMEN = (
    DIRECTORIO_SALIDAS
    / "tablas"
    / "tabla_02_resumen_cohorte.csv"
)
VALIDACION_GRUPOS = (
    DIRECTORIO_SALIDAS
    / "resultados"
    / "control_calidad"
    / "validacion_grupos_HER2_por_linea.csv"
)


# Valores obtenidos en la ejecución
VALORES_ESPERADOS = {
    "lineas": 22,
    "compuestos": 1383,
    "observaciones": 17253,
    "p20": 4.592383670461582,
    "p80": 7.386692396713945,
    "grupos": {"HER2_low": 5, "HER2_mid": 12, "HER2_high": 5},
}

COLUMNAS_IDENTIFICACION = [
    "depmap_id",
    "cell_line_display_name",
    "lineage_1",
    "lineage_2",
    "lineage_3",
]


def comprobar_archivo(ruta: Path) -> None:
    """Detiene la ejecución si falta un archivo de entrada."""
    if not ruta.is_file():
        raise FileNotFoundError(f"No se encontró el archivo requerido:\n{ruta}")


def comprobar_columnas(
    datos: pd.DataFrame,
    columnas: list[str],
    nombre_archivo: str,
) -> None:
    """Comprueba que una tabla contenga todas las columnas necesarias."""
    faltantes = sorted(set(columnas).difference(datos.columns))
    if faltantes:
        raise ValueError(
            f"En {nombre_archivo} faltan estas columnas: {', '.join(faltantes)}"
        )


def clasificar_grupo_her2(valor: float, p20: float, p80: float) -> str:
    """Clasifica una línea según los percentiles 20 y 80 de ERBB2."""
    if valor <= p20:
        return "HER2_low"
    if valor >= p80:
        return "HER2_high"
    return "HER2_mid"


def validar_resultado(
    datos: pd.DataFrame,
    expresion_por_linea: pd.DataFrame,
    p20: float,
    p80: float,
) -> None:
    """Verifica que el resultado coincida con la ejecución final del TFM."""
    observados = {
        "lineas": datos["depmap_id"].nunique(),
        "compuestos": datos["compound"].nunique(),
        "observaciones": len(datos),
    }

    errores = []
    for metrica in ("lineas", "compuestos", "observaciones"):
        if observados[metrica] != VALORES_ESPERADOS[metrica]:
            errores.append(
                f"{metrica}: {observados[metrica]} "
                f"(esperado: {VALORES_ESPERADOS[metrica]})"
            )

    if not isclose(p20, VALORES_ESPERADOS["p20"], abs_tol=1e-10):
        errores.append(f"P20: {p20} (esperado: {VALORES_ESPERADOS['p20']})")
    if not isclose(p80, VALORES_ESPERADOS["p80"], abs_tol=1e-10):
        errores.append(f"P80: {p80} (esperado: {VALORES_ESPERADOS['p80']})")

    grupos = expresion_por_linea["HER2_group"].value_counts().to_dict()
    if grupos != VALORES_ESPERADOS["grupos"]:
        errores.append(
            f"distribución de grupos: {grupos} "
            f"(esperada: {VALORES_ESPERADOS['grupos']})"
        )

    duplicados = datos.duplicated(["depmap_id", "compound"]).sum()
    if duplicados:
        errores.append(f"pares línea–compuesto duplicados: {duplicados}")

    if errores:
        raise RuntimeError(
            "Los datos no coinciden con la ejecución definitiva del TFM:\n- "
            + "\n- ".join(errores)
        )


def main() -> None:
    """Ejecuta la integración de los datos transcriptómicos y farmacológicos."""
    comprobar_archivo(ARCHIVO_AUC)
    comprobar_archivo(ARCHIVO_EXPRESION)

    datos_auc = pd.read_csv(ARCHIVO_AUC, low_memory=False)
    datos_expresion = pd.read_csv(ARCHIVO_EXPRESION, low_memory=False)

    comprobar_columnas(
        datos_auc,
        COLUMNAS_IDENTIFICACION,
        ARCHIVO_AUC.name,
    )
    comprobar_columnas(
        datos_expresion,
        COLUMNAS_IDENTIFICACION + ["ERBB2"],
        ARCHIVO_EXPRESION.name,
    )

    print("Dimensiones iniciales:")
    print(f"- Sensibilidad farmacológica: {datos_auc.shape}")
    print(f"- Expresión génica: {datos_expresion.shape}")

    # Seleccionar exclusivamente modelos de cáncer de mama.
    datos_auc = datos_auc.loc[datos_auc["lineage_1"].eq("Breast")].copy()
    datos_expresion = datos_expresion.loc[
        datos_expresion["lineage_1"].eq("Breast")
    ].copy()

    # Transformar PRISM de formato ancho a formato largo.
    auc_largo = datos_auc.melt(
        id_vars=COLUMNAS_IDENTIFICACION,
        var_name="compound",
        value_name="AUC",
    )
    auc_largo["AUC"] = pd.to_numeric(auc_largo["AUC"], errors="coerce")
    auc_largo = auc_largo.dropna(subset=["depmap_id", "compound", "AUC"])

    # Conservar una sola medición de ERBB2 por línea celular.
    expresion_erbb2 = datos_expresion[
        COLUMNAS_IDENTIFICACION + ["ERBB2"]
    ].copy()
    expresion_erbb2 = expresion_erbb2.rename(columns={"ERBB2": "ERBB2_expr"})
    expresion_erbb2["ERBB2_expr"] = pd.to_numeric(
        expresion_erbb2["ERBB2_expr"], errors="coerce"
    )
    expresion_erbb2 = expresion_erbb2.dropna(
        subset=["depmap_id", "ERBB2_expr"]
    )

    valores_por_linea = expresion_erbb2.groupby("depmap_id")["ERBB2_expr"].nunique()
    lineas_conflictivas = valores_por_linea[valores_por_linea > 1]
    if not lineas_conflictivas.empty:
        raise ValueError(
            "Existen líneas con más de un valor de ERBB2: "
            + ", ".join(lineas_conflictivas.index.astype(str))
        )

    expresion_erbb2 = expresion_erbb2.drop_duplicates(subset="depmap_id")

    # Integrar expresión de ERBB2 y AUC mediante el identificador DepMap.
    dataset_maestro = auc_largo.merge(
        expresion_erbb2[["depmap_id", "ERBB2_expr"]],
        on="depmap_id",
        how="inner",
        validate="many_to_one",
    )
    dataset_maestro = dataset_maestro.dropna(subset=["AUC", "ERBB2_expr"])
    if not dataset_maestro[["AUC", "ERBB2_expr"]].map(
        lambda x: float("-inf") < x < float("inf")
    ).all().all():
        raise ValueError("La integración contiene valores no finitos.")

    pares_duplicados = dataset_maestro.duplicated(
        ["depmap_id", "compound"], keep=False
    )
    if pares_duplicados.any():
        raise ValueError(
            "Se encontraron pares línea–compuesto duplicados después de la integración."
        )

    # Los percentiles se calculan sobre líneas únicas, no sobre observaciones.
    expresion_por_linea = dataset_maestro[
        ["depmap_id", "ERBB2_expr"]
    ].drop_duplicates(subset="depmap_id")
    p20 = expresion_por_linea["ERBB2_expr"].quantile(0.20)
    p80 = expresion_por_linea["ERBB2_expr"].quantile(0.80)

    dataset_maestro["HER2_group"] = dataset_maestro["ERBB2_expr"].apply(
        clasificar_grupo_her2,
        args=(p20, p80),
    )
    dataset_maestro = dataset_maestro.sort_values(
        ["depmap_id", "compound"]
    ).reset_index(drop=True)

    validacion_grupos = (
        dataset_maestro[
            COLUMNAS_IDENTIFICACION + ["ERBB2_expr", "HER2_group"]
        ]
        .drop_duplicates(subset="depmap_id")
        .sort_values("ERBB2_expr")
        .reset_index(drop=True)
    )

    validar_resultado(dataset_maestro, validacion_grupos, p20, p80)

    # Crear las carpetas necesarias y guardar las salidas.
    for ruta in (ARCHIVO_MAESTRO, TABLA_RESUMEN, VALIDACION_GRUPOS):
        ruta.parent.mkdir(parents=True, exist_ok=True)

    dataset_maestro.to_csv(ARCHIVO_MAESTRO, index=False)

    tabla_resumen = pd.DataFrame(
        {
            "Característica": [
                "Líneas celulares de cáncer de mama incluidas",
                "Compuestos evaluados",
                "Observaciones totales tras integración",
                "Percentil 20 de expresión de ERBB2",
                "Percentil 80 de expresión de ERBB2",
                "Líneas celulares HER2-low",
                "Líneas celulares HER2-mid",
                "Líneas celulares HER2-high",
                "Observaciones HER2-low",
                "Observaciones HER2-mid",
                "Observaciones HER2-high",
            ],
            "Valor": [
                dataset_maestro["depmap_id"].nunique(),
                dataset_maestro["compound"].nunique(),
                len(dataset_maestro),
                round(p20, 3),
                round(p80, 3),
                (validacion_grupos["HER2_group"] == "HER2_low").sum(),
                (validacion_grupos["HER2_group"] == "HER2_mid").sum(),
                (validacion_grupos["HER2_group"] == "HER2_high").sum(),
                (dataset_maestro["HER2_group"] == "HER2_low").sum(),
                (dataset_maestro["HER2_group"] == "HER2_mid").sum(),
                (dataset_maestro["HER2_group"] == "HER2_high").sum(),
            ],
        }
    )
    tabla_resumen.to_csv(TABLA_RESUMEN, index=False)
    validacion_grupos.to_csv(VALIDACION_GRUPOS, index=False)

    print("\nValidación correcta:")
    print(f"- Líneas celulares: {dataset_maestro['depmap_id'].nunique()}")
    print(f"- Compuestos: {dataset_maestro['compound'].nunique()}")
    print(f"- Observaciones: {len(dataset_maestro)}")
    print(f"- Percentil 20: {p20:.12f}")
    print(f"- Percentil 80: {p80:.12f}")
    print("- Grupos:")
    print(validacion_grupos["HER2_group"].value_counts().to_string())
    print(f"\nDataset maestro guardado en: {ARCHIVO_MAESTRO}")


if __name__ == "__main__":
    main()
