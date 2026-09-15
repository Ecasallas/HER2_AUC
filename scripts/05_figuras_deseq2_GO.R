#!/usr/bin/env Rscript
# Figuras 08–11 y tabla 06 a partir de las salidas verificadas del script 03.
# No vuelve a ejecutar DESeq2 ni enrichGO. No fija un número de términos GO.
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
})
args_script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(args_script)) stop("Ejecute este archivo mediante Rscript.")
script_path <- normalizePath(sub("^--file=", "", args_script[[1]]), mustWork = TRUE)
base_dir <- normalizePath(Sys.getenv("TFM_DIR", dirname(dirname(script_path))), mustWork = TRUE)
salidas_dir <- file.path(base_dir, "comprobacion_reproducibilidad")
dge_dir <- file.path(salidas_dir, "resultados", "expresion_diferencial")
go_dir <- file.path(salidas_dir, "resultados", "enriquecimiento_funcional")
apoyo_dir <- file.path(salidas_dir, "resultados", "tablas_apoyo")
control_dir <- file.path(salidas_dir, "resultados", "control_calidad")
fig_dir <- file.path(salidas_dir, "figuras")
tablas_dir <- file.path(salidas_dir, "tablas")
for (folder in c(fig_dir, tablas_dir, apoyo_dir)) dir.create(folder, recursive = TRUE, showWarnings = FALSE)

# La procedencia evita reutilizar los CSV anteriores con otro universo.
provenance_file <- file.path(control_dir, "procedencia_DESeq2_GO.rds")
if (!file.exists(provenance_file)) {
  stop("Ejecute primero 03_deseq2_y_enriquecimiento_GO.R: falta la procedencia verificable.")
}
prov <- readRDS(provenance_file)
paths <- file.path(salidas_dir, prov$archivos)
if (length(paths) != 5L || any(!file.exists(paths)) ||
    !identical(unname(tools::md5sum(paths)), prov$md5)) {
  stop("Las salidas cambiaron o están incompletas. Vuelva a ejecutar el script 03.")
}
master_file <- file.path(salidas_dir, "datos_procesados", "master_her2_auc.csv")
if (!file.exists(master_file) ||
    unname(tools::md5sum(master_file)) != prov$entradas$md5[[1]]) {
  stop("El maestro cambió desde la ejecución del script 03.")
}

deseq_file <- file.path(dge_dir, "deseq2_HER2_high_vs_low_resultados_prefiltrado.csv")
go_files <- c(BP = file.path(go_dir, "GO_proceso_biologico.csv"),
              MF = file.path(go_dir, "GO_funcion_molecular.csv"),
              CC = file.path(go_dir, "GO_componente_celular.csv"))
res <- read_csv(deseq_file, show_col_types = FALSE)
if (!all(c("gene", "log2FoldChange", "padj") %in% names(res))) stop("Columnas DESeq2 incompletas.")
if (any(res$padj < 0 | res$padj > 1, na.rm = TRUE)) stop("FDR fuera de [0,1].")
res <- res %>% filter(is.finite(log2FoldChange), !is.na(padj)) %>%
  mutate(Simbolo = str_remove(gene, "\\s*\\([0-9]+\\)$"),
         menos_log10_FDR = -log10(pmax(padj, .Machine$double.xmin)),
         Grupo = case_when(padj >= 0.05 ~ "No significativo",
                           log2FoldChange > 0 ~ "Mayor expresión en HER2-high",
                           log2FoldChange < 0 ~ "Mayor expresión en HER2-low",
                           TRUE ~ "Cambio nulo"))
n_sig <- sum(res$padj < 0.05)
n_pos <- sum(res$padj < 0.05 & res$log2FoldChange > 0)
n_neg <- sum(res$padj < 0.05 & res$log2FoldChange < 0)
if (!nrow(res)) stop("No hay genes representables.")

save_plot <- function(p, name, width = 10, height = 7) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = width, height = height,
         units = "in", dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = width, height = height,
         units = "in", bg = "white")
}
labels <- res %>% filter(Simbolo %in% c("ERBB2", "FGF2", "BCAS1", "TFF1", "TFF3", "ZEB1", "COL5A2", "ANPEP"))
volcano <- ggplot(res, aes(log2FoldChange, menos_log10_FDR, color = Grupo)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey45") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey65") +
  geom_label_repel(data = labels, aes(label = Simbolo), color = "black", size = 3,
                   seed = 20260915, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = c("No significativo" = "#BDBDBD",
                               "Mayor expresión en HER2-high" = "#16A6B6",
                               "Mayor expresión en HER2-low" = "#A85AC7",
                               "Cambio nulo" = "#777777")) +
  labs(title = "Expresión diferencial: HER2-high frente a HER2-low",
       subtitle = sprintf("%s genes con FDR < 0,05: %s positivos y %s negativos", n_sig, n_pos, n_neg),
       x = "Cambio de expresión (log2FoldChange)", y = expression(-log[10](FDR)),
       color = NULL, caption = "Grupos transcriptómicos. Significación definida por FDR < 0,05, sin umbral adicional de cambio.") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.text = element_text(size = 9),
        plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
save_plot(volcano, "figura_08_volcano_ERBB2", 11, 7.5)

go_terms_es <- c(
  # Biological Process
  "extracellular structure organization" = "Organización de estructuras extracelulares",
  "extracellular matrix organization" = "Organización de la matriz extracelular",
  "external encapsulating structure organization" = "Organización de estructuras envolventes externas",
  "connective tissue development" = "Desarrollo del tejido conectivo",
  "gliogenesis" = "Gliogénesis",
  "ear development" = "Desarrollo del oído",
  "myeloid leukocyte activation" = "Activación de leucocitos mieloides",
  "inner ear development" = "Desarrollo del oído interno",
  "cellular response to vascular endothelial growth factor stimulus" = "Respuesta celular al estímulo del factor de crecimiento endotelial vascular",
  "negative regulation of neuron projection development" = "Regulación negativa del desarrollo de proyecciones neuronales",
  "wound healing" = "Cicatrización",
  "regulation of body fluid levels" = "Regulación de los niveles de líquidos corporales",
  "negative regulation of cell projection organization" = "Regulación negativa de la organización de proyecciones celulares",
  "neuromuscular junction development" = "Desarrollo de la unión neuromuscular",
  "response to peptide hormone" = "Respuesta a hormonas peptídicas",
  
  # Molecular Function
  "extracellular matrix structural constituent" = "Componente estructural de la matriz extracelular",
  "growth factor binding" = "Unión a factores de crecimiento",
  "integrin binding" = "Unión a integrinas",
  "ligand-gated sodium channel activity" = "Actividad de canales de sodio regulados por ligando",
  "platelet-derived growth factor binding" = "Unión al factor de crecimiento derivado de plaquetas",
  "vascular endothelial growth factor receptor activity" = "Actividad del receptor del factor de crecimiento endotelial vascular",
  "ephrin receptor activity" = "Actividad del receptor de efrinas",
  "transmembrane-ephrin receptor activity" = "Actividad transmembrana del receptor de efrinas",
  "collagen receptor activity" = "Actividad del receptor de colágeno",
  "transmembrane receptor protein tyrosine kinase activity" = "Actividad tirosina quinasa de receptores transmembrana",
  "neurexin family protein binding" = "Unión a proteínas de la familia de las neurexinas",
  "GTPase activity" = "Actividad GTPasa",
  "GPI-linked ephrin receptor activity" = "Actividad del receptor de efrinas unido a GPI",
  "hepatocyte growth factor receptor activity" = "Actividad del receptor del factor de crecimiento hepatocitario",
  "insulin receptor activity" = "Actividad del receptor de insulina",
  
  # Cellular Component
  "collagen-containing extracellular matrix" = "Matriz extracelular con colágeno",
  "apical plasma membrane" = "Membrana plasmática apical",
  "transport vesicle" = "Vesícula de transporte",
  "apical part of cell" = "Región apical de la célula",
  "adherens junction" = "Unión adherente",
  "focal adhesion" = "Adhesión focal",
  "exocytic vesicle" = "Vesícula exocítica",
  "cell-substrate junction" = "Unión célula-sustrato",
  "endoplasmic reticulum lumen" = "Lumen del retículo endoplásmico",
  "apical junction complex" = "Complejo de unión apical",
  "transport vesicle membrane" = "Membrana de vesícula de transporte",
  "collagen trimer" = "Trímero de colágeno",
  "fibrillar collagen trimer" = "Trímero de colágeno fibrilar",
  "banded collagen fibril" = "Fibrilla de colágeno con bandas",
  "synaptic vesicle" = "Vesícula sináptica"
)
go_terms_es <- c(go_terms_es,
  "cell-cell junction" = "Unión célula-célula",
  "blood circulation" = "Circulación sanguínea",
  "endothelial cell migration" = "Migración de células endoteliales",
  "regulation of smooth muscle cell proliferation" = "Regulación de la proliferación de células musculares lisas")

# Traducciones por ID: se utiliza el diccionario del repositorio si está presente.
translation_file <- file.path(base_dir, "resultados", "tablas_apoyo", "traducciones_terminos_GO_figuras.csv")
translations <- data.frame(ID = character(), Description_es = character())
if (file.exists(translation_file)) {
  translations <- read_csv(translation_file, show_col_types = FALSE)
  if (!all(c("ID", "Description_es") %in% names(translations))) {
    stop("El diccionario debe contener ID y Description_es.")
  }
  if (anyDuplicated(translations$ID)) stop("Hay IDs GO duplicados en el diccionario.")
}

ontology_names <- c(BP = "Procesos biológicos", MF = "Funciones moleculares", CC = "Componentes celulares")
figure_names <- c(BP = "figura_09_GO_BP", MF = "figura_10_GO_MF", CC = "figura_11_GO_CC")
representative_tables <- list()
pending_translations <- list()
used_translations <- list()

parse_ratio <- function(x) {
  parts <- str_split_fixed(x, "/", 2)
  numerator <- as.numeric(parts[, 1]); denominator <- as.numeric(parts[, 2])
  if (anyNA(numerator) || anyNA(denominator) || any(denominator <= 0) ||
      any(numerator < 0 | numerator > denominator)) stop("GeneRatio/BgRatio no válido.")
  numerator / denominator
}

seleccion_tabla6 <- list(
  BP = c("extracellular structure organization", "extracellular matrix organization", "cellular response to vascular endothelial growth factor stimulus"),
  MF = c("extracellular matrix structural constituent", "integrin binding", "growth factor binding"),
  CC = c("collagen-containing extracellular matrix", "cell-cell junction", "transport vesicle")
)
for (ont in names(go_files)) {
  go <- read_csv(go_files[[ont]], show_col_types = FALSE)
  needed <- c("ID", "Description", "GeneRatio", "BgRatio", "p.adjust", "qvalue", "Count")
  if (!all(needed %in% names(go))) stop("Columnas GO incompletas para ", ont)
  if (nrow(go)) {
    if (anyNA(go$p.adjust) || anyNA(go$qvalue) || any(go$p.adjust < 0 | go$p.adjust >= 0.05) || any(go$qvalue >= 0.05 | go$qvalue < 0)) {
      stop("La tabla GO debe contener exclusivamente términos con FDR < 0,05: ", ont)
    }
    bg <- as.integer(sub(".*/", "", go$BgRatio))
    if (anyNA(bg) || any(bg > prov$n_universo)) stop("Universo GO incompatible: ", ont)
  }
  if (!nrow(go)) {
    p <- ggplot() + annotate("text", x = 0, y = 0,
                             label = "No se obtuvieron términos con FDR < 0,05", size = 5) +
      labs(title = paste("Enriquecimiento GO:", ontology_names[[ont]])) + theme_void()
    save_plot(p, figure_names[[ont]])
    representative_tables[[ont]] <- data.frame(Ontologia = character(), ID = character(),
      Termino_GO = character(), GeneRatio = character(), FDR = numeric())
    next
  }
  top <- go %>% arrange(p.adjust, ID) %>% slice_head(n = 15) %>%
    left_join(translations %>% select(ID, Description_es), by = "ID") %>%
    mutate(Termino = coalesce(na_if(Description_es, ""), unname(go_terms_es[Description]), Description),
           Ratio = parse_ratio(GeneRatio))
  pending_translations[[ont]] <- top %>% filter(Termino == Description) %>%
    select(ID, Description) %>% mutate(Description_es = "")
  used_translations[[ont]] <- top %>% select(ID, Description, Termino) %>% mutate(Ontologia = ont)
  # Selección de términos de la tabla 6 depositada.
  representative_tables[[ont]] <- go %>% filter(Description %in% seleccion_tabla6[[ont]]) %>%
    left_join(translations %>% select(ID, Description_es), by = "ID") %>%
    mutate(Termino = coalesce(na_if(Description_es, ""), unname(go_terms_es[Description]), Description)) %>%
    arrange(match(Description, seleccion_tabla6[[ont]])) %>%
    transmute(Ontologia = ont, ID, Termino_GO = Termino, GeneRatio, FDR = p.adjust)
  if (nrow(representative_tables[[ont]]) != 3L) stop("Falta algún término de la tabla 6 para ", ont)
  top <- top %>% arrange(Ratio, ID) %>%
    mutate(Etiqueta = paste0(str_wrap(Termino, 48), " [", ID, "]"),
           Etiqueta = factor(Etiqueta, levels = unique(Etiqueta)))
  p <- ggplot(top, aes(Ratio, Etiqueta, size = Count, color = p.adjust)) +
    geom_point(alpha = 0.95) +
    scale_color_gradient(low = "#B2182B", high = "#2166AC", name = "FDR") +
    scale_size_continuous(name = "Número de genes", range = c(3, 9)) +
    labs(title = paste("Enriquecimiento GO:", ontology_names[[ont]]),
         subtitle = sprintf("%d términos de menor FDR de un total de %d significativos", nrow(top), nrow(go)),
         x = "Proporción de genes (GeneRatio)", y = NULL,
         caption = "Universo: genes que superaron el prefiltrado de DESeq2. Genes diferenciales en ambas direcciones.") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 9), plot.margin = margin(10, 15, 10, 10))
  save_plot(p, figure_names[[ont]], 12, max(5, nrow(top) * 0.5 + 2))
}
write_csv(bind_rows(representative_tables), file.path(tablas_dir, "tabla_06_terminos_representativos_GO.csv"))
pending <- bind_rows(pending_translations) %>% distinct()
if (!ncol(pending)) pending <- data.frame(ID = character(), Description = character(), Description_es = character())
write_csv(pending, file.path(apoyo_dir, "traducciones_GO_pendientes.csv"))
write_csv(bind_rows(used_translations), file.path(apoyo_dir, "terminos_GO_representados.csv"))
if (nrow(pending)) warning("Algunos términos mantienen su nombre oficial en inglés. Complete traducciones_terminos_GO_figuras.csv y repita el script 05.")
capture.output(sessionInfo(), file = file.path(control_dir, "sessionInfo_figuras_R.txt"))
cat("Figuras 08–11 y tabla 06 generadas desde resultados con procedencia verificada.\n")
