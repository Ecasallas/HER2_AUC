#!/usr/bin/env Rscript

# ============================================================================
# TFM HER2-AUC: expresion diferencial con DESeq2 y enriquecimiento GO
# Universo GO corregido: genes que superaron el prefiltrado de DESeq2
# ============================================================================
#
# Este script reproduce y documenta la parte transcriptomica del TFM:
#   1. Recupera del master definitivo las 5 lineas de expresion baja de ERBB2
#      y las 5 lineas de expresion alta de ERBB2.
#   2. Selecciona sus conteos crudos de RNA-seq.
#   3. Ejecuta DESeq2 con el contraste alta frente a baja.
#   4. Selecciona genes con FDR < 0,05.
#   5. Ejecuta el analisis de sobrerrepresentacion de Gene Ontology (GO)
#      por separado para BP, MF y CC.
#
# IMPORTANTE SOBRE EL UNIVERSO DE GO
# ----------------------------------
# enrichGO() utiliza explicitamente como universo los genes que superaron el
# prefiltrado de DESeq2 (conteo >= 10 en al menos 2 de las 10 lineas). Los
# identificadores ENTREZ se extraen directamente de los nombres de las columnas
# de la matriz de DepMap, que siguen el formato "SIMBOLO (ENTREZID)".
#
# Entradas esperadas (relativas a TFM_DIR):
#   data_processed/master_her2_auc.csv
#   data_raw/OmicsExpressionRawReadCountHumanProteinCodingGenes.csv
#
# Salidas:
#   results_dge/reproducible_final_universe_prefiltrado/
#
# Se utiliza una carpeta nueva para no sobrescribir el analisis GO anterior.
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# ----------------------------------------------------------------------------
# 0. Configuracion
# ----------------------------------------------------------------------------
DEFAULT_TFM_DIR <-
  "/Users/estefanialejandracasallassamper/Desktop/TFM_UNIR/TFM_HER2_AUC"

BASE_DIR <- Sys.getenv("TFM_DIR", unset = DEFAULT_TFM_DIR)

MASTER_FILE <- file.path(
  BASE_DIR,
  "data_processed",
  "master_her2_auc.csv"
)

COUNTS_FILE <- file.path(
  BASE_DIR,
  "data_raw",
  "OmicsExpressionRawReadCountHumanProteinCodingGenes.csv"
)

OUT_DIR <- file.path(
  BASE_DIR,
  "results_dge",
  "reproducible_final_universe_prefiltrado"
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Estas expectativas corresponden a resultados que no deben cambiar al
# corregir el universo de GO. Los numeros de terminos GO no se fijan porque
# pueden cambiar al utilizar el fondo de genes realmente evaluados.
EXPECTED <- list(
  n_master_lines = 22L,
  n_low = 5L,
  n_high = 5L,
  p20 = 4.592383670461582,
  p80 = 7.386692396713945,
  n_genes_after_prefilter = 15774L,
  n_significant = 1708L,
  n_positive = 893L,
  n_negative = 815L
)

MIN_COUNT <- 10L
MIN_SAMPLES <- 2L
FDR_CUTOFF <- 0.05

# Si TRUE, una diferencia en los resultados fundamentales detiene el script
# para evitar mezclar accidentalmente versiones de datos.
STRICT_VALIDATION <- TRUE

# ----------------------------------------------------------------------------
# 1. Funciones auxiliares
# ----------------------------------------------------------------------------
require_file <- function(path) {
  if (!file.exists(path)) {
    stop("No se encontro el archivo requerido:\n", path, call. = FALSE)
  }
}

check_expected <- function(observed, expected, label, tolerance = NULL) {
  matches <- if (is.null(tolerance)) {
    identical(as.integer(observed), as.integer(expected))
  } else {
    isTRUE(all.equal(
      as.numeric(observed),
      as.numeric(expected),
      tolerance = tolerance
    ))
  }

  if (!matches) {
    message_text <- paste0(
      label, ": valor observado = ", observed,
      "; valor esperado = ", expected, "."
    )
    if (STRICT_VALIDATION) {
      stop(message_text, call. = FALSE)
    } else {
      warning(message_text, call. = FALSE)
    }
  }
}

standardize_id_column <- function(df, object_name) {
  candidates <- c("depmap_id", "ModelID", "model_id", "...1")
  found <- candidates[candidates %in% names(df)]

  if (length(found) == 0L) {
    stop(
      "No se encontro un identificador de modelo en ", object_name,
      ". Se esperaba una de estas columnas: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }

  if (found[[1]] != "depmap_id") {
    names(df)[names(df) == found[[1]]] <- "depmap_id"
  }
  df
}

filter_default_entries <- function(df) {
  default_candidates <- c(
    "IsDefaultEntryForModel",
    "is_default_entry",
    "IsDefaultForModel"
  )
  default_col <- default_candidates[default_candidates %in% names(df)]

  if (length(default_col) == 0L) {
    message(
      "La matriz de conteos no contiene IsDefaultEntryForModel; ",
      "se comprobaran directamente los duplicados por depmap_id."
    )
    return(df)
  }

  default_col <- default_col[[1]]
  normalized <- toupper(trimws(as.character(df[[default_col]])))
  keep <- normalized %in% c("YES", "TRUE", "1", "Y", "SI", "SÍ")

  if (!any(keep)) {
    stop(
      "Existe la columna ", default_col,
      ", pero no se reconocio ninguna entrada predeterminada.",
      call. = FALSE
    )
  }

  message(
    "Se conservaron ", sum(keep),
    " entradas con ", default_col, " = Yes."
  )
  df[keep, , drop = FALSE]
}

# Extrae el identificador numerico situado entre los ultimos parentesis.
# Ejemplo: "FGF2 (2247)" -> "2247".
extract_entrez_id <- function(x) {
  extracted <- sub(
    "^.*\\(([0-9]+)\\)\\s*$",
    "\\1",
    as.character(x)
  )
  valid <- grepl("^[0-9]+$", extracted)
  extracted[!valid] <- NA_character_
  extracted
}

extract_gene_symbol <- function(x) {
  trimws(sub("\\s*\\([0-9]+\\)\\s*$", "", as.character(x)))
}

extract_bg_denominator <- function(go_df) {
  if (nrow(go_df) == 0L || !("BgRatio" %in% names(go_df))) {
    return(NA_integer_)
  }
  as.integer(sub(".*/", "", go_df$BgRatio[[1]]))
}

extract_gene_denominator <- function(go_df) {
  if (nrow(go_df) == 0L || !("GeneRatio" %in% names(go_df))) {
    return(NA_integer_)
  }
  as.integer(sub(".*/", "", go_df$GeneRatio[[1]]))
}

# ----------------------------------------------------------------------------
# 2. Carga y validacion del master definitivo
# ----------------------------------------------------------------------------
require_file(MASTER_FILE)
require_file(COUNTS_FILE)

master <- read_csv(MASTER_FILE, show_col_types = FALSE)
master <- standardize_id_column(master, "master_her2_auc.csv")

required_master <- c("depmap_id", "ERBB2_expr", "HER2_group")
missing_master <- setdiff(required_master, names(master))
if (length(missing_master) > 0L) {
  stop(
    "Faltan columnas en el master: ",
    paste(missing_master, collapse = ", "),
    call. = FALSE
  )
}

master <- master %>%
  mutate(ERBB2_expr = as.numeric(ERBB2_expr)) %>%
  filter(!is.na(depmap_id), !is.na(ERBB2_expr), !is.na(HER2_group))

line_table <- master %>%
  dplyr::select(depmap_id, ERBB2_expr, HER2_group) %>%
  distinct()

conflicting_lines <- line_table %>%
  count(depmap_id) %>%
  filter(n != 1L)

if (nrow(conflicting_lines) > 0L) {
  stop(
    "Existen lineas con mas de un valor de ERBB2 o mas de un grupo en el master.",
    call. = FALSE
  )
}

check_expected(
  n_distinct(line_table$depmap_id),
  EXPECTED$n_master_lines,
  "Numero de lineas del master"
)

p20 <- quantile(line_table$ERBB2_expr, 0.20, names = FALSE)
p80 <- quantile(line_table$ERBB2_expr, 0.80, names = FALSE)

check_expected(p20, EXPECTED$p20, "Percentil 20", tolerance = 1e-10)
check_expected(p80, EXPECTED$p80, "Percentil 80", tolerance = 1e-10)

meta <- line_table %>%
  filter(HER2_group %in% c("HER2_low", "HER2_high")) %>%
  mutate(
    HER2_group = factor(
      HER2_group,
      levels = c("HER2_low", "HER2_high")
    )
  ) %>%
  arrange(HER2_group, ERBB2_expr)

group_counts <- table(meta$HER2_group)
check_expected(
  group_counts[["HER2_low"]],
  EXPECTED$n_low,
  "Lineas de expresion baja"
)
check_expected(
  group_counts[["HER2_high"]],
  EXPECTED$n_high,
  "Lineas de expresion alta"
)

write_csv(
  meta %>% mutate(HER2_group = as.character(HER2_group)),
  file.path(OUT_DIR, "lineas_incluidas_DESeq2.csv")
)

# ----------------------------------------------------------------------------
# 3. Preparacion de la matriz de conteos crudos
# ----------------------------------------------------------------------------
counts <- read_csv(COUNTS_FILE, show_col_types = FALSE)
counts <- standardize_id_column(counts, "matriz de conteos")
counts <- filter_default_entries(counts)

counts_selected <- counts %>%
  filter(depmap_id %in% meta$depmap_id)

missing_lines <- setdiff(meta$depmap_id, counts_selected$depmap_id)
if (length(missing_lines) > 0L) {
  stop(
    "Faltan en la matriz de conteos estas lineas seleccionadas:\n",
    paste(missing_lines, collapse = "\n"),
    call. = FALSE
  )
}

duplicated_ids <- counts_selected$depmap_id[duplicated(counts_selected$depmap_id)]
if (length(duplicated_ids) > 0L) {
  stop(
    "Quedo mas de una entrada de conteos para: ",
    paste(unique(duplicated_ids), collapse = ", "),
    ". No se seleccionara una fila de forma arbitraria.",
    call. = FALSE
  )
}

counts_selected <- counts_selected[
  match(meta$depmap_id, counts_selected$depmap_id),
  ,
  drop = FALSE
]

# En la matriz de DepMap, las columnas genicas siguen el patron
# "SIMBOLO (ENTREZID)". Asi se evita interpretar metadatos como genes.
gene_cols <- names(counts_selected)[
  str_detect(names(counts_selected), "\\s\\([0-9]+\\)$")
]

if (length(gene_cols) == 0L) {
  stop(
    "No se detectaron columnas genicas con el patron 'SIMBOLO (ENTREZID)'.",
    call. = FALSE
  )
}

if (anyDuplicated(gene_cols)) {
  stop("La matriz contiene nombres de genes duplicados.", call. = FALSE)
}

numeric_counts <- counts_selected %>%
  dplyr::select(all_of(gene_cols)) %>%
  mutate(across(everything(), as.numeric))

count_matrix_samples_by_genes <- as.matrix(numeric_counts)

if (anyNA(count_matrix_samples_by_genes)) {
  stop(
    "La conversion numerica produjo valores ausentes en los conteos.",
    call. = FALSE
  )
}

if (any(!is.finite(count_matrix_samples_by_genes))) {
  stop("La matriz contiene valores no finitos.", call. = FALSE)
}

if (any(count_matrix_samples_by_genes < 0)) {
  stop("La matriz contiene conteos negativos.", call. = FALSE)
}

if (any(abs(
  count_matrix_samples_by_genes - round(count_matrix_samples_by_genes)
) > 1e-8)) {
  stop(
    "La matriz contiene valores no enteros. No se redondearan ni truncaran ",
    "automaticamente porque DESeq2 requiere conteos crudos o estimados ",
    "apropiados.",
    call. = FALSE
  )
}

count_matrix <- t(round(count_matrix_samples_by_genes))
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- gene_cols
colnames(count_matrix) <- counts_selected$depmap_id

meta2 <- meta[match(colnames(count_matrix), meta$depmap_id), , drop = FALSE]
meta2 <- as.data.frame(meta2)
rownames(meta2) <- meta2$depmap_id

if (!identical(colnames(count_matrix), rownames(meta2))) {
  stop("El orden de la matriz y los metadatos no coincide.", call. = FALSE)
}

# ----------------------------------------------------------------------------
# 4. Expresion diferencial con DESeq2
# ----------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = meta2,
  design = ~HER2_group
)

# Prefiltrado: conteo >= 10 en al menos 2 de las 10 lineas celulares.
keep_genes <- rowSums(counts(dds) >= MIN_COUNT) >= MIN_SAMPLES
dds <- dds[keep_genes, ]

check_expected(
  nrow(dds),
  EXPECTED$n_genes_after_prefilter,
  "Genes conservados despues del prefiltrado"
)

dds <- DESeq(dds)

# HER2_low es la referencia. Un log2 fold change positivo indica mayor
# expresion en el grupo de expresion alta de ERBB2.
res <- results(
  dds,
  contrast = c("HER2_group", "HER2_high", "HER2_low"),
  alpha = FDR_CUTOFF,
  independentFiltering = TRUE
)

# No se aplica shrinkage a los log2 fold changes presentados.
res_tbl <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj, pvalue)

sig_tbl <- res_tbl %>%
  filter(!is.na(padj), padj < FDR_CUTOFF)

n_positive <- sum(sig_tbl$log2FoldChange > 0, na.rm = TRUE)
n_negative <- sum(sig_tbl$log2FoldChange < 0, na.rm = TRUE)

check_expected(
  nrow(sig_tbl),
  EXPECTED$n_significant,
  "Genes con FDR < 0,05"
)
check_expected(
  n_positive,
  EXPECTED$n_positive,
  "Genes con cambio positivo"
)
check_expected(
  n_negative,
  EXPECTED$n_negative,
  "Genes con cambio negativo"
)

write_csv(res_tbl, file.path(OUT_DIR, "DESeq2_resultados_completos.csv"))
write_csv(sig_tbl, file.path(OUT_DIR, "DESeq2_genes_FDR_menor_0_05.csv"))
write_csv(head(res_tbl, 20), file.path(OUT_DIR, "DESeq2_top20_por_FDR.csv"))

# ----------------------------------------------------------------------------
# 5. Identificadores para Gene Ontology
# ----------------------------------------------------------------------------
# Los ENTREZID se extraen directamente del formato "SIMBOLO (ENTREZID)".
gene_map <- sig_tbl %>%
  transmute(
    gene_original = gene,
    SYMBOL = extract_gene_symbol(gene),
    ENTREZID = extract_entrez_id(gene)
  )

genes_without_entrez <- gene_map %>%
  filter(is.na(ENTREZID) | ENTREZID == "")

if (nrow(genes_without_entrez) > 0L) {
  stop(
    "No se pudo extraer un identificador ENTREZ de ",
    nrow(genes_without_entrez),
    " genes significativos. Ejemplos: ",
    paste(head(genes_without_entrez$gene_original, 5L), collapse = ", "),
    call. = FALSE
  )
}

gene_map <- gene_map %>%
  filter(SYMBOL != "", ENTREZID != "") %>%
  distinct(ENTREZID, .keep_all = TRUE)

if (nrow(gene_map) != EXPECTED$n_significant) {
  stop(
    "Se obtuvieron ", nrow(gene_map),
    " identificadores ENTREZ unicos; se esperaban ",
    EXPECTED$n_significant,
    ". Revise si existen identificadores duplicados.",
    call. = FALSE
  )
}

entrez_significant <- unique(gene_map$ENTREZID)

write_csv(
  gene_map,
  file.path(OUT_DIR, "genes_significativos_convertidos_ENTREZID.csv")
)

# Universo: todos los genes que superaron el prefiltrado de DESeq2, no solo
# los genes diferenciales.
universe_map <- tibble(gene_original = rownames(dds)) %>%
  transmute(
    gene_original = gene_original,
    SYMBOL = extract_gene_symbol(gene_original),
    ENTREZID = extract_entrez_id(gene_original)
  )

universe_without_entrez <- universe_map %>%
  filter(is.na(ENTREZID) | ENTREZID == "")

if (nrow(universe_without_entrez) > 0L) {
  stop(
    "No se pudo extraer un identificador ENTREZ de ",
    nrow(universe_without_entrez),
    " genes del universo. Ejemplos: ",
    paste(head(universe_without_entrez$gene_original, 5L), collapse = ", "),
    call. = FALSE
  )
}

universe_map <- universe_map %>%
  filter(SYMBOL != "", ENTREZID != "") %>%
  distinct(ENTREZID, .keep_all = TRUE)

entrez_universe <- unique(universe_map$ENTREZID)

if (length(entrez_universe) == 0L) {
  stop("El universo GO quedo vacio.", call. = FALSE)
}

missing_from_universe <- setdiff(entrez_significant, entrez_universe)
if (length(missing_from_universe) > 0L) {
  stop(
    "Hay ", length(missing_from_universe),
    " genes significativos que no estan presentes en el universo GO.",
    call. = FALSE
  )
}

write_csv(
  universe_map,
  file.path(OUT_DIR, "universo_GO_genes_prefiltrados.csv")
)

message("Genes significativos para GO: ", length(entrez_significant))
message("Genes unicos del universo GO: ", length(entrez_universe))

# ----------------------------------------------------------------------------
# 6. Enriquecimiento GO con universo explicito
# ----------------------------------------------------------------------------
run_enrich_go <- function(ontology) {
  enrichGO(
    gene = entrez_significant,
    universe = entrez_universe,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = FDR_CUTOFF,
    qvalueCutoff = FDR_CUTOFF,
    readable = TRUE
  )
}

go_bp <- run_enrich_go("BP")
go_mf <- run_enrich_go("MF")
go_cc <- run_enrich_go("CC")

go_bp_df <- as.data.frame(go_bp)
go_mf_df <- as.data.frame(go_mf)
go_cc_df <- as.data.frame(go_cc)

# Los conteos GO pueden diferir de los 254/41/34 del analisis anterior porque
# ahora se utiliza el universo de genes prefiltrados. Se informan los valores
# observados sin imponer una validacion contra los numeros antiguos.
if (nrow(go_bp_df) == 0L) {
  warning("No se obtuvieron terminos GO significativos para BP.", call. = FALSE)
}
if (nrow(go_mf_df) == 0L) {
  warning("No se obtuvieron terminos GO significativos para MF.", call. = FALSE)
}
if (nrow(go_cc_df) == 0L) {
  warning("No se obtuvieron terminos GO significativos para CC.", call. = FALSE)
}

write_csv(
  go_bp_df,
  file.path(OUT_DIR, "GO_Biological_Process_enrichment.csv")
)
write_csv(
  go_mf_df,
  file.path(OUT_DIR, "GO_Molecular_Function_enrichment.csv")
)
write_csv(
  go_cc_df,
  file.path(OUT_DIR, "GO_Cellular_Component_enrichment.csv")
)

# ----------------------------------------------------------------------------
# 7. Resumen y trazabilidad
# ----------------------------------------------------------------------------
summary_table <- tibble(
  metrica = c(
    "Lineas del master",
    "Lineas de expresion baja",
    "Lineas de expresion alta",
    "Percentil 20",
    "Percentil 80",
    "Genes despues del prefiltrado",
    "Genes con FDR < 0,05",
    "Genes con cambio positivo",
    "Genes con cambio negativo",
    "Genes significativos con ENTREZID",
    "Genes unicos del universo GO",
    "Terminos GO BP",
    "Terminos GO MF",
    "Terminos GO CC",
    "Fondo efectivo GO BP",
    "Fondo efectivo GO MF",
    "Fondo efectivo GO CC",
    "Genes de entrada anotados en BP",
    "Genes de entrada anotados en MF",
    "Genes de entrada anotados en CC"
  ),
  valor = c(
    n_distinct(line_table$depmap_id),
    unname(group_counts[["HER2_low"]]),
    unname(group_counts[["HER2_high"]]),
    p20,
    p80,
    nrow(dds),
    nrow(sig_tbl),
    n_positive,
    n_negative,
    length(entrez_significant),
    length(entrez_universe),
    nrow(go_bp_df),
    nrow(go_mf_df),
    nrow(go_cc_df),
    extract_bg_denominator(go_bp_df),
    extract_bg_denominator(go_mf_df),
    extract_bg_denominator(go_cc_df),
    extract_gene_denominator(go_bp_df),
    extract_gene_denominator(go_mf_df),
    extract_gene_denominator(go_cc_df)
  )
)

write_csv(summary_table, file.path(OUT_DIR, "resumen_validacion.csv"))

capture.output(
  sessionInfo(),
  file = file.path(OUT_DIR, "sessionInfo.txt")
)

cat("\nAnalisis completado con universo GO explicito.\n")
cat("- Lineas: 5 de expresion baja y 5 de expresion alta.\n")
cat("- Genes tras prefiltrado:", nrow(dds), "\n")
cat("- Genes con FDR < 0,05:", nrow(sig_tbl), "\n")
cat("- Cambios positivos:", n_positive, "\n")
cat("- Cambios negativos:", n_negative, "\n")
cat("- Genes unicos del universo GO:", length(entrez_universe), "\n")
cat("- Terminos GO BP:", nrow(go_bp_df), "\n")
cat("- Terminos GO MF:", nrow(go_mf_df), "\n")
cat("- Terminos GO CC:", nrow(go_cc_df), "\n")
cat("- Salidas guardadas en:", OUT_DIR, "\n")
