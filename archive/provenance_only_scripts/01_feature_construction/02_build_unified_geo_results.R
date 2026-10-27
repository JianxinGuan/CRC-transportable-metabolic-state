invisible(Sys.setlocale("LC_CTYPE", "Chinese"))
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

out <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "geo_stage15")
dir.create(file.path(out, "unified_results"), recursive = TRUE, showWarnings = FALSE)

empty_col <- function(n) rep(NA_character_, n)

std_deseq <- function(file, dataset, contrast, cell_line, treatment_type,
                      platform, evidence_level, integration_role,
                      gene_map = NULL, gene_id_is_symbol = FALSE) {
  x <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  gene_id <- as.character(x[[1]])
  if (gene_id_is_symbol) {
    symbol <- gene_id
  } else if (!is.null(gene_map)) {
    symbol <- unname(gene_map[gene_id])
  } else {
    symbol <- empty_col(nrow(x))
  }
  data.frame(
    dataset = dataset, contrast = contrast, contrast_family = treatment_type,
    cell_line = cell_line, gene_id = gene_id, gene_symbol = symbol,
    entrez_id = empty_col(nrow(x)), log2FC = x$log2FoldChange,
    standard_error = x$lfcSE, test_statistic = x$stat,
    P_Value = x$pvalue, adj_P_Value = x$padj,
    analysis_method = "DESeq2", platform = platform,
    evidence_level = evidence_level, integration_role = integration_role,
    source_file = basename(file), stringsAsFactors = FALSE
  )
}

std_deseq_tsv <- function(file, dataset, contrast, cell_line, treatment_type,
                          platform, evidence_level, integration_role) {
  x <- read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
  data.frame(
    dataset = dataset, contrast = contrast, contrast_family = treatment_type,
    cell_line = cell_line, gene_id = x$gene_ensembl, gene_symbol = x$gene_symbol,
    entrez_id = empty_col(nrow(x)), log2FC = x$log2FoldChange,
    standard_error = x$lfcSE, test_statistic = x$stat,
    P_Value = x$pvalue, adj_P_Value = x$padj,
    analysis_method = "DESeq2", platform = platform,
    evidence_level = evidence_level, integration_role = integration_role,
    source_file = basename(file), stringsAsFactors = FALSE
  )
}

std_limma_csv <- function(file, dataset, contrast, cell_line, treatment_type,
                          platform, evidence_level, integration_role) {
  x <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  gene <- as.character(x[[1]])
  se <- abs(x$logFC / x$t)
  se[!is.finite(se)] <- NA_real_
  data.frame(
    dataset = dataset, contrast = contrast, contrast_family = treatment_type,
    cell_line = cell_line, gene_id = gene, gene_symbol = gene,
    entrez_id = empty_col(nrow(x)), log2FC = x$logFC,
    standard_error = se, test_statistic = x$t,
    P_Value = x$P.Value, adj_P_Value = x$adj.P.Val,
    analysis_method = "limma", platform = platform,
    evidence_level = evidence_level, integration_role = integration_role,
    source_file = basename(file), stringsAsFactors = FALSE
  )
}

std_limma_tsv <- function(file, dataset, contrast, cell_line, treatment_type,
                          platform, evidence_level, integration_role) {
  x <- read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
  se <- abs(x$logFC / x$t)
  se[!is.finite(se)] <- NA_real_
  ent <- if ("entrez_id" %in% colnames(x)) as.character(x$entrez_id) else empty_col(nrow(x))
  data.frame(
    dataset = dataset, contrast = contrast, contrast_family = treatment_type,
    cell_line = cell_line, gene_id = x$gene_symbol, gene_symbol = x$gene_symbol,
    entrez_id = ent, log2FC = x$logFC,
    standard_error = se, test_statistic = x$t,
    P_Value = x$P.Value, adj_P_Value = x$adj.P.Val,
    analysis_method = "limma", platform = platform,
    evidence_level = evidence_level, integration_role = integration_role,
    source_file = basename(file), stringsAsFactors = FALSE
  )
}

res <- list()
add <- function(x) res[[length(res) + 1L]] <<- x

# RNA-seq cohorts.
add(std_deseq(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "s09/GSE200546_DESeq2_all_results.csv"), "GSE200546",
              "pH6.6_vs_pH7.4", "CLS1", "acute_acidosis",
              "RNA-seq expected counts", "primary", "primary_dataset_level",
              gene_id_is_symbol = TRUE))

g116 <- read.delim(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "s08/expression_matrix/GSE116035_raw_counts.tsv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
map116 <- setNames(g116$gene_symbol, g116$gene_ensembl)
add(std_deseq(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "s11/overall_pH_effect.csv"), "GSE116035",
              "pH6.5_vs_pH7.4_adjusted_cell_line", "FaDu;HCT-116;SiHa",
              "acute_acidosis", "RNA-seq processed counts", "primary",
              "primary_dataset_level", gene_map = map116))
for (cl in c("FaDu", "HCT-116", "SiHa")) {
  safe <- gsub("[^A-Za-z0-9]+", "_", cl)
  add(std_deseq_tsv(file.path(out, "GSE116035_stratified_DESeq2",
                              paste0("GSE116035_", safe, "_pH6.5_vs_pH7.4_all.tsv")),
                    "GSE116035", paste0(cl, "_pH6.5_vs_pH7.4"), cl,
                    "acute_acidosis", "RNA-seq processed counts", "sensitivity",
                    "cell_line_sensitivity_only"))
}

ens_keys <- unique(sub("\\..*$", "", as.character(read.csv(
  file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "s12/PANC1_DESeq2_results.csv"), check.names = FALSE)[[1]])))
ens_symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = ens_keys, column = "SYMBOL",
                                    keytype = "ENSEMBL", multiVals = "first")
for (cl in c("MDA-MB-231", "PANC1", "MCF7")) {
  fn <- if (cl == "MDA-MB-231") "MDA-MB-231_DESeq2_results.csv" else
    paste0(cl, "_DESeq2_results.csv")
  add(std_deseq(file.path(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "s12"), fn), "GSE152345", "pH6.5_vs_pH7.6", cl,
                "acute_acidosis", "RNA-seq gene counts", "primary",
                "primary_cell_line_context", gene_map = ens_symbol))
}

# GPL571 gene-level rerun.
for (cn in c("LA_vs_Ctrl", "Hypoxia_vs_Ctrl", "LA_Hypoxia_vs_Ctrl")) {
  family <- if (cn == "LA_vs_Ctrl") "lactic_acidosis" else if (cn == "Hypoxia_vs_Ctrl") "hypoxia" else "lactic_acidosis_hypoxia"
  role <- if (cn == "LA_vs_Ctrl") "primary_dataset_level" else "auxiliary_stress_evidence"
  level <- if (cn == "LA_vs_Ctrl") "primary" else "auxiliary"
  add(std_limma_tsv(file.path(out, "GSE29406_GPL571_gene_level",
                              paste0("GSE29406_", cn, "_gene_level_all.tsv")),
                    "GSE29406", cn, "MCF7", family, "GPL571", level, role))
}

# GSE300768 uses a documented backup annotation and remains secondary.
add(std_limma_csv(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "s13/GSE300768_Acidosis_vs_Control_all.csv"), "GSE300768",
                  "pH6.4_vs_pH7.4", "PANC1", "chronic_acidosis",
                  "GPL23126 backup annotation", "secondary",
                  "secondary_pending_official_annotation"))

# GSE9649 planned contrasts; experiment block was included in the source design.
g964 <- data.frame(
  code = c("LA_vs_Ctrl", "Acidosis_vs_Ctrl", "Hypoxia_vs_Ctrl", "LA_Hypoxia_vs_Ctrl", "Lactosis_vs_Ctrl"),
  family = c("lactic_acidosis", "acute_acidosis", "hypoxia", "lactic_acidosis_hypoxia", "lactosis"),
  level = c("primary", "primary", "auxiliary", "auxiliary", "auxiliary"), stringsAsFactors = FALSE)
for (i in seq_len(nrow(g964))) {
  role <- if (g964$level[i] == "primary") "primary_dataset_level" else "auxiliary_stress_evidence"
  add(std_limma_csv(file.path(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "s13"), paste0("GSE9649_", g964$code[i], "_all.csv")),
                    "GSE9649", g964$code[i], "HMEC", g964$family[i], "GPL570",
                    g964$level[i], role))
}

# GSE300765: two acidosis contrasts are primary; hypoxia is auxiliary.
g300765 <- data.frame(
  code = c("chronic_acidosis_vs_chronic_neutral", "acute_acidosis_vs_acute_neutral", "hypoxia_vs_normoxia"),
  family = c("chronic_acidosis", "acute_acidosis", "hypoxia"),
  level = c("primary", "primary", "auxiliary"), stringsAsFactors = FALSE)
for (i in seq_len(nrow(g300765))) {
  role <- if (g300765$level[i] == "primary") "primary_within_dataset_contrast" else "auxiliary_stress_evidence"
  add(std_limma_tsv(file.path(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "s14/limma_results"), paste0(g300765$code[i], "_all.tsv")),
                    "GSE300765", g300765$code[i], "U87MG", g300765$family[i],
                    "GPL10558", g300765$level[i], role))
}

unified <- do.call(rbind, res)
rownames(unified) <- NULL

# Complete missing Entrez IDs from stable human symbols, preserving source IDs.
symbols <- unique(unified$gene_symbol[!is.na(unified$gene_symbol) & nzchar(unified$gene_symbol)])
symbol_entrez <- suppressMessages(AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))
need_entrez <- is.na(unified$entrez_id) | !nzchar(unified$entrez_id)
unified$entrez_id[need_entrez] <- unname(symbol_entrez[unified$gene_symbol[need_entrez]])
unified$direction <- ifelse(is.na(unified$log2FC), "not_estimable",
                            ifelse(unified$log2FC > 0, "up_in_treatment",
                                   ifelse(unified$log2FC < 0, "down_in_treatment", "no_change")))
unified$pass_abs_log2FC2_P005 <- !is.na(unified$log2FC) & !is.na(unified$P_Value) &
  abs(unified$log2FC) >= 2 & unified$P_Value < 0.05
unified$pass_abs_log2FC2_FDR005 <- !is.na(unified$log2FC) & !is.na(unified$adj_P_Value) &
  abs(unified$log2FC) >= 2 & unified$adj_P_Value < 0.05
unified$gene_symbol_valid <- !is.na(unified$gene_symbol) & nzchar(unified$gene_symbol)

con <- gzfile(file.path(out, "unified_results", "GEO_all_contrasts_unified.tsv.gz"), "wt")
write.table(unified, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
close(con)
sig <- unified[unified$pass_abs_log2FC2_P005 & unified$gene_symbol_valid, ]
write.table(sig, file.path(out, "unified_results", "GEO_nominal_logFC2_P005_unified.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

split_key <- interaction(unified$dataset, unified$contrast, unified$cell_line,
                         drop = TRUE, lex.order = TRUE)
summ <- do.call(rbind, lapply(split(unified, split_key), function(z) data.frame(
  dataset = z$dataset[1], contrast = z$contrast[1], contrast_family = z$contrast_family[1],
  cell_line = z$cell_line[1], analysis_method = z$analysis_method[1], platform = z$platform[1],
  evidence_level = z$evidence_level[1], integration_role = z$integration_role[1],
  tested_rows = nrow(z), valid_gene_symbols = sum(z$gene_symbol_valid),
  up_logFC2_P005 = sum(z$log2FC >= 2 & z$P_Value < 0.05, na.rm = TRUE),
  down_logFC2_P005 = sum(z$log2FC <= -2 & z$P_Value < 0.05, na.rm = TRUE),
  up_logFC2_FDR005 = sum(z$log2FC >= 2 & z$adj_P_Value < 0.05, na.rm = TRUE),
  down_logFC2_FDR005 = sum(z$log2FC <= -2 & z$adj_P_Value < 0.05, na.rm = TRUE)
)))
rownames(summ) <- NULL
write.table(summ, file.path(out, "unified_results", "GEO_contrast_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

audit <- data.frame(
  dataset = c("GSE200546", "GSE116035", "GSE152345", "GSE29406", "GSE300768",
              "GSE9649", "GSE300765", "GSE215085", "GSE315655"),
  assay = c("RNA-seq", "RNA-seq", "RNA-seq", "microarray", "microarray",
            "microarray", "microarray", "processed fold-change", "unknown/incomplete"),
  latest_status = c("DESeq2 complete", "overall and 3 stratified DESeq2 complete",
                    "3 cell-line DESeq2 complete", "GPL571 annotated gene-level limma complete",
                    "limma complete with backup annotation", "GPL570 gene-level limma complete",
                    "GPL10558 gene-level limma complete", "fold-change tables only",
                    "series matrix incomplete"),
  primary_analysis = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE),
  integration_use = c("primary", "primary overall; strata sensitivity", "primary by cell line",
                      "primary LA; hypoxia/combined auxiliary", "secondary sensitivity only",
                      "primary acidosis/LA; other stresses auxiliary",
                      "primary chronic/acute; hypoxia auxiliary", "auxiliary only", "excluded"),
  exclusion_or_caution = c(
    "RSEM expected counts rounded for DESeq2; 24H/0H auxiliary samples excluded",
    "processed counts rounded; cell-line heterogeneity explicitly assessed",
    "do not pool three cell lines as biological replicates",
    "duplicate probes collapsed by maximum IQR after hgu133a2.db annotation",
    "GPL23126 annotation is a documented GitHub backup, not official platform export",
    "multiple treatments share study context; experiment block retained",
    "actual cell line is U87MG, correcting obsolete audit entry that stated PANC1",
    "no sample-level expression matrix; cannot enter DESeq2/limma meta-analysis",
    "3 KB metadata/incomplete matrix; no valid expression analysis"),
  stringsAsFactors = FALSE
)
write.table(audit, file.path(out, "final_cohort_eligibility_audit.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

writeLines(c(
  "Stage 15 final audit and unified GEO differential-expression package",
  "Positive log2FC always means higher expression in the named treatment/stress group versus its control.",
  "Primary integration must not count overall and cell-line sensitivity contrasts from GSE116035 as independent evidence.",
  "Multiple contrasts within GSE9649 and GSE300765 are not independent cohorts.",
  "GSE300768 remains secondary until official GPL23126 annotation sensitivity validation is available.",
  "GSE215085 is auxiliary only; GSE315655 is excluded.",
  "The full unified table is gzip-compressed; the nominal abs(log2FC)>=2, P<0.05 subset is uncompressed.",
  paste0("org.Hs.eg.db=", as.character(packageVersion("org.Hs.eg.db")))
), file.path(out, "README_STAGE15.txt"))

writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo_unification.txt"))
cat("Unified rows:", nrow(unified), " significant rows:", nrow(sig), " contrasts:", nrow(summ), "\n")
