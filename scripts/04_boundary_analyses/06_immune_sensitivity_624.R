args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("Usage: script patient.tsv mcp.tsv xcell.tsv c7.tsv.gz c7_diff.tsv checkpoint_expression.tsv output_dir")
}

patient_file <- args[1]
mcp_file <- args[2]
xcell_file <- args[3]
c7_file <- args[4]
c7_old_diff_file <- args[5]
checkpoint_file <- args[6]
out <- args[7]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_literal_matrix <- function(f) {
  con <- if (grepl("[.]gz$", f, ignore.case = TRUE)) gzfile(f) else f
  zz <- readLines(con, warn = FALSE)
  sp <- strsplit(zz, "`t", fixed = TRUE)
  ids <- sp[[1]]
  vals <- lapply(sp[-1], function(v) as.numeric(v[-1]))
  m <- do.call(rbind, vals)
  rownames(m) <- vapply(sp[-1], `[`, "", 1)
  colnames(m) <- ids
  m
}

wilcox_matrix <- function(mat, patients, delta_name = "delta_C2_minus_C1") {
  ids <- intersect(colnames(mat), patients$sample_id)
  cl <- patients$cluster[match(ids, patients$sample_id)]
  ans <- do.call(rbind, lapply(seq_len(nrow(mat)), function(i) {
    y <- as.numeric(mat[i, ids])
    ok <- is.finite(y) & !is.na(cl)
    y <- y[ok]
    g <- cl[ok]
    if (length(unique(g)) < 2L) {
      return(data.frame(feature = rownames(mat)[i], n = length(y), C1_n = sum(g == 1),
        C2_n = sum(g == 2), C1_mean = if (any(g == 1)) mean(y[g == 1]) else NA_real_,
        C2_mean = if (any(g == 2)) mean(y[g == 2]) else NA_real_, delta = NA_real_,
        p_value = NA_real_, stringsAsFactors = FALSE))
    }
    w <- wilcox.test(y ~ factor(g, levels = c(1, 2)), exact = FALSE)
    data.frame(feature = rownames(mat)[i], n = length(y), C1_n = sum(g == 1),
      C2_n = sum(g == 2), C1_mean = mean(y[g == 1]), C2_mean = mean(y[g == 2]),
      delta = mean(y[g == 2]) - mean(y[g == 1]), p_value = w$p.value,
      stringsAsFactors = FALSE)
  }))
  ans$FDR <- p.adjust(ans$p_value, method = "BH")
  names(ans)[names(ans) == "delta"] <- delta_name
  ans
}

patients <- read.delim(patient_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(patients) == 624L, length(unique(patients$patient_barcode_fixed)) == 624L)
patients$cluster <- as.integer(patients$cluster)

mcp <- read_literal_matrix(mcp_file)
xcell <- read_literal_matrix(xcell_file)
c7 <- read_literal_matrix(c7_file)

mcp_res <- wilcox_matrix(mcp, patients)
xcell_res <- wilcox_matrix(xcell, patients)
c7_res <- wilcox_matrix(c7, patients)
write.table(mcp_res, file.path(out, "MCPcounter_differential_624.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(xcell_res, file.path(out, "xCell_differential_624.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(c7_res, file.path(out, "C7_differential_624.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Frozen manual mapping used in the original cross-method comparison.
mapping <- data.frame(
  MCPcounter_cell = c("CD8 T cells", "B lineage", "Endothelial cells", "Fibroblasts", "T cells", "NK cells", "Neutrophils"),
  xCell_cell = c("CD8+ T-cells", "pro B-cells", "ly Endothelial cells", "Fibroblasts", "Th1 cells", "NK cells", "Neutrophils"),
  stringsAsFactors = FALSE
)
cross <- merge(mapping, transform(mcp_res, MCPcounter_cell = feature), by = "MCPcounter_cell", all.x = TRUE)
cross <- merge(cross, transform(xcell_res, xCell_cell = feature), by = "xCell_cell", all.x = TRUE,
               suffixes = c("_MCP", "_xCell"))
cross$direction_concordant <- sign(cross$delta_C2_minus_C1_MCP) == sign(cross$delta_C2_minus_C1_xCell)
cross$at_least_one_FDR05 <- cross$FDR_MCP < 0.05 | cross$FDR_xCell < 0.05
cross$robust_rule <- cross$direction_concordant & cross$at_least_one_FDR05
write.table(cross, file.path(out, "MCP_xCell_cross_validation_624.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Freeze theme membership from the original 200-signature result by keyword rules used for reporting.
theme_patterns <- list(
  CD8_T = "CD8", T_cell = "TCELL|T_CELL|T-CELL|T CELL|CD4|TREG", B_cell = "BCELL|B_CELL", NK = "NKT|NK_CELL|NKCELL",
  Myeloid = "MONOCYTE|MACROPHAGE|MYELOID|NEUTROPHIL|DENDRITIC", Interferon = "INTERFERON|IFN",
  Immune_suppression = "TREG|REGULATORY|SUPPRESS|CHECKPOINT|PD1|PD_1|CTLA4"
)
theme <- do.call(rbind, lapply(names(theme_patterns), function(nm) {
  q <- c7_res[grepl(theme_patterns[[nm]], c7_res$feature, ignore.case = TRUE), , drop = FALSE]
  data.frame(theme = nm, c7_n = nrow(q), c7_up_in_C2 = sum(q$delta_C2_minus_C1 > 0),
    c7_down_in_C2 = sum(q$delta_C2_minus_C1 < 0), c7_fdr05 = sum(q$FDR < 0.05),
    c7_median_delta_C2_minus_C1 = if (nrow(q)) median(q$delta_C2_minus_C1) else NA_real_)
}))
write.table(theme, file.path(out, "C7_theme_summary_624.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Checkpoint/effector expression table is gene x sample on the original TPM scale used in stage 37.
chk <- read.delim(checkpoint_file, row.names = 1, check.names = FALSE)
chk_res <- wilcox_matrix(as.matrix(chk), patients)
write.table(chk_res, file.path(out, "checkpoint_effector_differential_624.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

old_lines <- strsplit(readLines(c7_old_diff_file, warn = FALSE), "`t", fixed = TRUE)
old_c7 <- as.data.frame(do.call(rbind, old_lines[-1]), stringsAsFactors = FALSE)
names(old_c7) <- c("feature", "old_delta", "old_p_value", "old_FDR")
old_c7[-1] <- lapply(old_c7[-1], as.numeric)
audit_c7 <- merge(old_c7, c7_res[, c("feature", "delta_C2_minus_C1", "p_value", "FDR")], by = "feature", all = TRUE)
names(audit_c7)[5:7] <- c("new_delta", "new_p_value", "new_FDR")
audit_c7$direction_changed <- sign(-audit_c7$old_delta) != sign(audit_c7$new_delta)
audit_c7$significance_changed <- (audit_c7$old_FDR < 0.05) != (audit_c7$new_FDR < 0.05)
write.table(audit_c7, file.path(out, "old_vs_new_C7_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

summary <- data.frame(
  analysis = c("MCP-counter", "xCell", "C7", "checkpoint_effector", "MCP_xCell_robust"),
  features = c(nrow(mcp_res), nrow(xcell_res), nrow(c7_res), nrow(chk_res), nrow(cross)),
  FDR05_or_robust = c(sum(mcp_res$FDR < 0.05), sum(xcell_res$FDR < 0.05),
    sum(c7_res$FDR < 0.05), sum(chk_res$FDR < 0.05), sum(cross$robust_rule)),
  stringsAsFactors = FALSE
)
write.table(summary, file.path(out, "immune_sensitivity_summary_624.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
capture.output(sessionInfo(), file = file.path(out, "R_sessionInfo.txt"))
