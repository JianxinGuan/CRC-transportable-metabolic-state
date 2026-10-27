invisible(Sys.setlocale("LC_CTYPE", "Chinese"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(AnnotationDbi)
  library(hgu133a2.db)
})

out <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "geo_stage15")
dir.create(file.path(out, "GSE116035_stratified_DESeq2"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "GSE29406_GPL571_gene_level"), recursive = TRUE, showWarnings = FALSE)

# GSE116035: independent within-cell-line models. The previous overall object
# cannot be safely subset because its design retains a one-level cell_line term.
counts_df <- read.delim(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "s08/expression_matrix/GSE116035_raw_counts.tsv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.delim(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "s08/sample_metadata/GSE116035_sample_metadata.tsv"),
                   stringsAsFactors = FALSE)
sample_cols <- meta$sample_id
stopifnot(all(sample_cols %in% colnames(counts_df)))
counts <- round(as.matrix(counts_df[, sample_cols, drop = FALSE]))
storage.mode(counts) <- "integer"
rownames(counts) <- counts_df$gene_ensembl

strat_summary <- list()
for (cl in unique(meta$cell_line)) {
  m <- meta[meta$cell_line == cl, , drop = FALSE]
  m <- m[match(colnames(counts[, m$sample_id, drop = FALSE]), m$sample_id), , drop = FALSE]
  rownames(m) <- m$sample_id
  m$pH <- relevel(factor(m$pH), ref = "pH7.4")
  cmat <- counts[, m$sample_id, drop = FALSE]
  dds <- DESeqDataSetFromMatrix(cmat, m, design = ~pH)
  dds <- dds[rowSums(counts(dds) > 0) >= 2, ]
  dds <- DESeq(dds, fitType = "parametric", sfType = "poscounts", quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("pH", "pH6.5", "pH7.4")))
  res$gene_ensembl <- rownames(res)
  symbol_map <- setNames(counts_df$gene_symbol, counts_df$gene_ensembl)
  res$gene_symbol <- unname(symbol_map[res$gene_ensembl])
  res <- res[, c("gene_ensembl", "gene_symbol", "baseMean", "log2FoldChange",
                 "lfcSE", "stat", "pvalue", "padj")]
  safe_cl <- gsub("[^A-Za-z0-9]+", "_", cl)
  write.table(res, file.path(out, "GSE116035_stratified_DESeq2",
                             paste0("GSE116035_", safe_cl, "_pH6.5_vs_pH7.4_all.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  strat_summary[[cl]] <- data.frame(
    dataset = "GSE116035", cell_line = cl, n_samples = nrow(m), tested_genes = nrow(res),
    nominal_LFC2_P005 = sum(abs(res$log2FoldChange) >= 2 & res$pvalue < 0.05, na.rm = TRUE),
    FDR_LFC2_005 = sum(abs(res$log2FoldChange) >= 2 & res$padj < 0.05, na.rm = TRUE)
  )
}
write.table(do.call(rbind, strat_summary),
            file.path(out, "GSE116035_stratified_DESeq2", "GSE116035_stratified_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# GSE29406: annotate GPL571 (Affymetrix HG-U133A 2.0) using the pinned
# hgu133a2.db package, then retain the most variable probe per gene symbol.
g294 <- read.delim(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "s10/GSE29406_expression_matrix.tsv"), check.names = FALSE,
                   stringsAsFactors = FALSE)
gmeta <- read.delim(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "s10/GSE29406_sample_metadata.tsv"), stringsAsFactors = FALSE)
probes <- g294$probe_id
raw_ann <- AnnotationDbi::select(hgu133a2.db, keys = unique(probes),
                                 columns = c("SYMBOL", "ENTREZID"), keytype = "PROBEID")
raw_ann <- raw_ann[order(raw_ann$PROBEID, is.na(raw_ann$SYMBOL), raw_ann$SYMBOL), ]
raw_ann$mapping_rank <- ave(seq_len(nrow(raw_ann)), raw_ann$PROBEID, FUN = seq_along)
raw_ann$mapping_multiplicity <- ave(raw_ann$PROBEID, raw_ann$PROBEID, FUN = length)
write.table(raw_ann, file.path(out, "GSE29406_GPL571_gene_level",
                              "GPL571_hgu133a2db_all_probe_mappings.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
ann <- raw_ann[raw_ann$mapping_rank == 1, ]
ann <- ann[match(probes, ann$PROBEID), ]

expr <- as.matrix(g294[, gmeta$sample_id, drop = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- probes
probe_iqr <- apply(expr, 1, IQR, na.rm = TRUE)
probe_ann <- data.frame(probe_id = probes, gene_symbol = ann$SYMBOL,
                        entrez_id = ann$ENTREZID,
                        mapping_multiplicity = ann$mapping_multiplicity,
                        probe_IQR = probe_iqr, stringsAsFactors = FALSE)
valid <- !is.na(probe_ann$gene_symbol) & nzchar(probe_ann$gene_symbol)
ord <- order(probe_ann$gene_symbol, -probe_ann$probe_IQR, probe_ann$probe_id)
selected <- ord[valid[ord] & !duplicated(probe_ann$gene_symbol[ord])]
probe_ann$selected_for_gene_level <- FALSE
probe_ann$selected_for_gene_level[selected] <- TRUE
write.table(probe_ann, file.path(out, "GSE29406_GPL571_gene_level",
                                "GSE29406_GPL571_probe_annotation_audit.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

gene_expr <- expr[selected, , drop = FALSE]
rownames(gene_expr) <- probe_ann$gene_symbol[selected]
gene_map <- probe_ann[selected, c("gene_symbol", "probe_id", "entrez_id", "probe_IQR")]
write.table(data.frame(gene_symbol = rownames(gene_expr), gene_expr, check.names = FALSE),
            file.path(out, "GSE29406_GPL571_gene_level", "GSE29406_gene_expression_matrix.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(gene_map, file.path(out, "GSE29406_GPL571_gene_level",
                               "GSE29406_selected_probe_per_gene.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

gmeta$group <- factor(gmeta$group)
design <- model.matrix(~0 + group, gmeta)
colnames(design) <- sub("group", "", colnames(design))
rownames(design) <- gmeta$sample_id
cts <- makeContrasts(
  LA_vs_Ctrl = lactic_acidosis - control,
  Hypoxia_vs_Ctrl = hypoxia - control,
  LA_Hypoxia_vs_Ctrl = lactic_acidosis_hypoxia - control,
  levels = design
)
fit <- eBayes(contrasts.fit(lmFit(gene_expr, design), cts), trend = TRUE, robust = TRUE)
limma_summary <- list()
for (cn in colnames(cts)) {
  tt <- topTable(fit, coef = cn, number = Inf, sort.by = "P")
  tt$gene_symbol <- rownames(tt)
  tt$entrez_id <- gene_map$entrez_id[match(tt$gene_symbol, gene_map$gene_symbol)]
  tt$selected_probe <- gene_map$probe_id[match(tt$gene_symbol, gene_map$gene_symbol)]
  tt <- tt[, c("gene_symbol", "entrez_id", "selected_probe", "logFC", "AveExpr",
               "t", "P.Value", "adj.P.Val", "B")]
  write.table(tt, file.path(out, "GSE29406_GPL571_gene_level",
                           paste0("GSE29406_", cn, "_gene_level_all.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  sig <- tt[abs(tt$logFC) >= 2 & tt$P.Value < 0.05, ]
  write.table(sig, file.path(out, "GSE29406_GPL571_gene_level",
                            paste0("GSE29406_", cn, "_gene_level_logFC2_P005.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  limma_summary[[cn]] <- data.frame(
    dataset = "GSE29406", contrast = cn, probes_input = nrow(expr),
    probes_annotated = sum(valid), unique_genes_tested = nrow(tt),
    nominal_LFC2_P005 = nrow(sig),
    FDR_LFC2_005 = sum(abs(tt$logFC) >= 2 & tt$adj.P.Val < 0.05, na.rm = TRUE)
  )
}
write.table(do.call(rbind, limma_summary),
            file.path(out, "GSE29406_GPL571_gene_level", "GSE29406_gene_level_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

writeLines(c(
  paste0("R=", R.version.string),
  paste0("DESeq2=", as.character(packageVersion("DESeq2"))),
  paste0("limma=", as.character(packageVersion("limma"))),
  paste0("AnnotationDbi=", as.character(packageVersion("AnnotationDbi"))),
  paste0("hgu133a2.db=", as.character(packageVersion("hgu133a2.db")))
), file.path(out, "analysis_package_versions.txt"))

cat("Stage 15 core analyses completed.\n")
