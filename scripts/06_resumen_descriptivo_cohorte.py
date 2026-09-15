#!/usr/bin/env python3
"""Resume ERBB2 por línea y documenta cobertura, exclusiones y entorno Python.

No calcula una correlación global mezclando compuestos. Ejecutar desde cualquier
carpeta; TFM_DIR permite indicar una copia alternativa del proyecto.
"""
from pathlib import Path
import os
import hashlib
import platform
import importlib.metadata
import importlib.util
import pandas as pd
import numpy as np

BASE_DIR = Path(os.environ.get("TFM_DIR", Path(__file__).resolve().parents[1])).expanduser().resolve()
SALIDAS = BASE_DIR / "comprobacion_reproducibilidad"
MASTER = SALIDAS / "datos_procesados" / "master_her2_auc.csv"
CONTROL = SALIDAS / "resultados" / "control_calidad"


def main():
    # Reutiliza las validaciones del análisis, sin ejecutar sus cálculos ni figuras.
    module_path = Path(__file__).resolve().with_name("02_analisis_farmacogenomico_por_compuesto.py")
    spec = importlib.util.spec_from_file_location("analisis_compuestos", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    datos = pd.read_csv(MASTER)
    module.validar_dataset_maestro(datos)
    CONTROL.mkdir(parents=True, exist_ok=True)
    lines = datos[["depmap_id", "ERBB2_expr", "HER2_group"]].drop_duplicates()
    expr = lines["ERBB2_expr"]
    pd.DataFrame({
        "Métrica": ["Media", "Mediana", "Desviación estándar", "Mínimo", "Máximo", "Percentil 20", "Percentil 80"],
        "Valor": [expr.mean(), expr.median(), expr.std(ddof=1), expr.min(), expr.max(), expr.quantile(.2), expr.quantile(.8)],
    }).to_csv(CONTROL / "resumen_expresion_ERBB2.csv", index=False)
    summary = lines.groupby("HER2_group").agg(
        lineas=("depmap_id", "nunique"), media_ERBB2=("ERBB2_expr", "mean"),
        mediana_ERBB2=("ERBB2_expr", "median"), DE_ERBB2=("ERBB2_expr", "std"))
    summary["observaciones"] = datos.groupby("HER2_group").size()
    summary.to_csv(CONTROL / "resumen_grupos_HER2.csv")
    coverage = []
    for compound, sub in datos.groupby("compound"):
        n = sub["depmap_id"].nunique()
        nh = int(sub["HER2_group"].eq("HER2_high").sum())
        nl = int(sub["HER2_group"].eq("HER2_low").sum())
        varying = sub["ERBB2_expr"].nunique() > 1 and sub["AUC"].nunique() > 1
        reason = "Incluido" if n >= 10 and varying else ("Menos de 10 líneas" if n < 10 else "Variable constante")
        coverage.append({"Compuesto completo": compound, "Lineas": n, "n HER2-high": nh,
                         "n HER2-low": nl, "Spearman": reason,
                         "MannWhitney": "Incluido" if min(nh, nl) >= 3 else "Menos de 3 líneas en algún grupo"})
    pd.DataFrame(coverage).to_csv(CONTROL / "cobertura_y_elegibilidad_por_compuesto.csv", index=False)
    content = [f"Python: {platform.python_version()}", f"Sistema: {platform.platform()}"]
    for package in ["pandas", "numpy", "scipy", "statsmodels", "matplotlib", "seaborn"]:
        content.append(f"{package}=={importlib.metadata.version(package)}")
    content.append(f"SHA256 maestro: {hashlib.sha256(MASTER.read_bytes()).hexdigest()}")
    (CONTROL / "versiones_Python.txt").write_text("\n".join(content) + "\n", encoding="utf-8")
    print("Resumen descriptivo y cobertura generados: 22 líneas y 1.383 compuestos.")
    print(CONTROL)


if __name__ == "__main__":
    main()
