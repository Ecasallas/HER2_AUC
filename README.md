# HER2_AUC

Análisis farmacogenómico de la expresión de ERBB2 y la respuesta a fármacos en líneas celulares de cáncer de mama mediante datos públicos de DepMap y PRISM.

Este repositorio contiene los datos utilizados, los scripts, los resultados, las tablas y las figuras generadas para el Trabajo Fin de Máster:

**Análisis farmacogenómico de la expresión de ERBB2 y la sensibilidad a fármacos en líneas celulares de cáncer de mama mediante datos de DepMap y PRISM**

**Autora:** Estefanía Alejandra Casallas Samper  
**Programa:** Máster Universitario en Bioinformática  
**Universidad:** Universidad Internacional de La Rioja (UNIR)  
**Año:** 2026

## Objetivo

El objetivo del estudio fue evaluar la asociación entre la expresión transcriptómica de ERBB2 y la sensibilidad farmacológica en líneas celulares de cáncer de mama.

Para ello, se integraron datos de expresión génica procedentes de DepMap/CCLE con datos de sensibilidad farmacológica del programa PRISM Repurposing.

La sensibilidad farmacológica se evaluó mediante el área bajo la curva dosis-respuesta (AUC). Un valor de AUC menor indica una mayor sensibilidad al compuesto, mientras que un valor mayor indica una menor sensibilidad.

## Diseño del estudio

Se realizó un estudio observacional, analítico y transversal basado en el análisis secundario de datos públicos.

La unidad de análisis farmacogenómico fue cada combinación única línea celular–compuesto. Cada observación incluyó la expresión de ERBB2 de la línea celular y el valor de AUC correspondiente al compuesto evaluado.

La cohorte integrada definitiva estuvo formada por:

- 22 líneas celulares de cáncer de mama.
- 1.383 compuestos.
- 17.253 observaciones línea–compuesto.
- 5 líneas HER2-low.
- 12 líneas HER2-mid.
- 5 líneas HER2-high.

Los grupos extremos de expresión de ERBB2 se definieron mediante los percentiles 20 y 80 de la expresión, medida en log2(TPM + 1):

- Percentil 20: 4,592.
- Percentil 80: 7,387.

Las denominaciones HER2-low, HER2-mid y HER2-high corresponden a grupos transcriptómicos definidos para este estudio y no equivalen a las categorías clínicas basadas en inmunohistoquímica o hibridación in situ.

## Análisis realizados

El flujo de trabajo incluyó los siguientes análisis:

1. Integración de los datos transcriptómicos y farmacológicos.
2. Control de calidad y validación de la cohorte.
3. Clasificación de las líneas celulares según la expresión de ERBB2.
4. Correlación de Spearman entre la expresión de ERBB2 y el AUC para cada compuesto.
5. Corrección por comparaciones múltiples mediante el procedimiento de Benjamini–Hochberg.
6. Comparación entre los grupos HER2-high y HER2-low mediante la prueba de Mann–Whitney.
7. Análisis de expresión diferencial con DESeq2.
8. Análisis de enriquecimiento funcional mediante Gene Ontology.
9. Generación reproducible de tablas y figuras.

## Resultados principales

La correlación de Spearman mostró asociaciones farmacológicas heterogéneas, cuya magnitud y dirección dependieron del compuesto analizado.

Dinaciclib fue el único compuesto que alcanzó el criterio principal de significación estadística, FDR inferior a 0,05, después del ajuste por comparaciones múltiples:

- Rho de Spearman: 0,778.
- Valor p: 3,30 × 10⁻⁵.
- FDR: 0,029.
- Número de líneas celulares evaluadas: 21.

La correlación positiva indica que una mayor expresión de ERBB2 se asoció con un AUC mayor y, por tanto, con una menor sensibilidad a dinaciclib.

También se identificaron señales exploratorias para inhibidores de HER/EGFR, como canertinib, poziotinib, lapatinib y tucatinib, además de alpelisib, un inhibidor de PI3Kα.

En la comparación entre los grupos HER2-high y HER2-low mediante la prueba de Mann–Whitney, ningún compuesto permaneció significativo después de la corrección por FDR. Por tanto, estos resultados deben interpretarse como tendencias exploratorias.

El análisis de expresión diferencial entre los grupos extremos incluyó:

- 10 líneas celulares: 5 HER2-high y 5 HER2-low.
- 15.774 genes después del prefiltrado.
- 1.708 genes diferencialmente expresados con FDR inferior a 0,05.
- 893 genes con `log2FoldChange` positivo.
- 815 genes con `log2FoldChange` negativo.

El contraste se definió como HER2-high frente a HER2-low. Los valores positivos de `log2FoldChange` indican mayor expresión en HER2-high y los negativos, menor expresión en este grupo.

ERBB2 presentó una mayor expresión en el grupo HER2-high, de acuerdo con la clasificación utilizada.

El enriquecimiento funcional mostró términos relacionados con la organización de la matriz extracelular, la unión a integrinas, la unión a factores de crecimiento, las uniones celulares y las vesículas de transporte.

## Estructura del repositorio

| Directorio o archivo | Contenido |
|---|---|
| `datos_originales/` | Archivos de entrada y subconjuntos de los datos originales obtenidos de DepMap y PRISM. |
| `datos_procesados/` | Dataset maestro generado después de integrar y filtrar los datos. |
| `scripts/` | Scripts utilizados para procesar los datos, ejecutar los análisis y generar las figuras. |
| `resultados/` | Resultados completos de los análisis farmacogenómicos, transcriptómicos y funcionales. |
| `tablas/` | Tablas principales utilizadas en el documento final y tablas complementarias. |
| `figuras/` | Figuras generadas para representar los resultados del estudio. |
| `LICENSE` | Licencia aplicable al código original del repositorio. |
| `README.md` | Descripción general del proyecto y guía del contenido del repositorio. |

## Datos originales

El directorio `datos_originales/` contiene los siguientes archivos:

- `Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_subsetted.csv`
- `Expression_Public_25Q3_subsetted.csv`
- `OmicsExpressionRawReadCountHumanProteinCodingGenes.csv.gz`

Los dos primeros archivos corresponden a subconjuntos de los datos originales, preparados para reproducir el análisis farmacogenómico de las líneas celulares de cáncer de mama.

El archivo `OmicsExpressionRawReadCountHumanProteinCodingGenes.csv.gz` contiene la matriz de conteos crudos de genes codificantes utilizada en el análisis de expresión diferencial con DESeq2. Se conserva comprimido en formato gzip para reducir su tamaño. Puede leerse directamente mediante `readr::read_csv()`, sin necesidad de descomprimirlo previamente.

El archivo principal obtenido durante la integración de los datos transcriptómicos y farmacológicos se encuentra en:

- `datos_procesados/master_her2_auc.csv`

Los datos transcriptómicos corresponden a DepMap Public 25Q3 y los datos farmacológicos proceden de PRISM Repurposing.

Los datos originales no son propiedad de la autora de este repositorio y no están cubiertos por su licencia MIT. Su uso, distribución y atribución se mantienen sujetos a las condiciones establecidas por sus proveedores originales.

## Organización de los resultados

El directorio `resultados/` está dividido en las siguientes categorías:

| Subdirectorio | Descripción |
|---|---|
| `control_calidad/` | Archivos de validación, resúmenes de la cohorte y controles de calidad. |
| `farmacogenomica/` | Resultados de las correlaciones de Spearman y comparaciones de Mann–Whitney. |
| `expresion_diferencial/` | Resultados completos y filtrados del análisis realizado con DESeq2. |
| `enriquecimiento_funcional/` | Resultados del análisis de enriquecimiento de Gene Ontology. |
| `tablas_apoyo/` | Archivos utilizados para generar tablas, figuras o resúmenes del estudio. |

Los resultados completos se mantienen separados de las tablas y figuras utilizadas en el documento final.

## Requisitos

El análisis se desarrolló principalmente con Python y R.

### Python

Principales bibliotecas utilizadas:

- pandas
- numpy
- scipy
- matplotlib
- seaborn

Las dependencias pueden instalarse con:

```bash
pip install pandas numpy scipy matplotlib seaborn
```

### R

**Versión utilizada:** R 4.5.0

Principales paquetes utilizados:

- DESeq2
- readr
- dplyr
- tibble
- stringr
- clusterProfiler
- org.Hs.eg.db

Los paquetes de Bioconductor pueden instalarse mediante:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(c(
    "DESeq2",
    "clusterProfiler",
    "org.Hs.eg.db"
))
```

Los demás paquetes pueden instalarse mediante:

```r
install.packages(c(
    "readr",
    "dplyr",
    "tibble",
    "stringr"
))
```

Estos comandos no fijan las versiones de los paquetes. Las diferencias entre versiones de software y bases de anotación pueden afectar a la reproducción exacta de los resultados, especialmente en el enriquecimiento funcional.

## Reproducción del análisis

Para reproducir el flujo de trabajo:

1. Descargar o clonar este repositorio.
2. Mantener la estructura original de directorios.
3. Instalar las dependencias de Python y R.
4. Comprobar que los archivos de entrada se encuentran en `datos_originales/`.
5. Comprobar que las rutas configuradas en los scripts corresponden a la estructura del repositorio y que el archivo de conteos se lee con su extensión `.csv.gz`.
6. Ejecutar los scripts de integración y control de calidad.
7. Ejecutar los análisis farmacogenómicos.
8. Ejecutar el análisis de expresión diferencial con DESeq2.
9. Ejecutar el enriquecimiento funcional.
10. Generar las tablas y figuras finales.

Para clonar el repositorio:

```bash
git clone https://github.com/Ecasallas/HER2_AUC.git
cd HER2_AUC
```

## Interpretación de los resultados

Este estudio tiene carácter exploratorio y se basa en modelos celulares in vitro.

El reducido número de líneas celulares de los grupos extremos limita la potencia estadística, especialmente en las comparaciones realizadas mediante la prueba de Mann–Whitney.

Los resultados representan asociaciones estadísticas y no permiten establecer relaciones causales ni aplicaciones clínicas directas. Su confirmación requiere el análisis de cohortes independientes y la validación en modelos experimentales adicionales.

## Fuentes de datos

Los datos utilizados proceden de los siguientes recursos:

- DepMap Portal y la versión DepMap Public 25Q3.
- PRISM Repurposing.
- Cancer Cell Line Encyclopedia (CCLE).

Estos recursos se integran en el contexto del proyecto Cancer Dependency Map del Broad Institute.

La utilización y redistribución de estos datos debe respetar las condiciones establecidas por sus proveedores originales.

## Licencia

El código original incluido en el directorio `scripts/` se distribuye bajo la licencia MIT:

```text
Copyright (c) 2026 Estefanía Alejandra Casallas Samper
```

La licencia MIT no se aplica a los datos de terceros contenidos en `datos_originales/` ni sustituye las condiciones de uso establecidas por DepMap, PRISM o el Broad Institute.

Los archivos incluidos en `datos_procesados/`, `resultados/`, `tablas/` y `figuras/` fueron generados a partir del análisis de datos públicos y se proporcionan con fines de transparencia y reproducibilidad académica. Su disponibilidad no elimina las condiciones que puedan resultar aplicables a los datos de origen.

Consulta el archivo `LICENSE` para obtener más información.

## Contacto

**Estefanía Alejandra Casallas Samper**

Repositorio: [github.com/Ecasallas/HER2_AUC](https://github.com/Ecasallas/HER2_AUC)
