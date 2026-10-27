invisible(Sys.setlocale("LC_CTYPE", "Chinese"))
suppressPackageStartupMessages(library(data.table))

out <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "geo_stage16")
dir.create(file.path(out, "01_feasibility"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "02_dataset_level_effects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "03_human_meta_results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "04_crispr_intersections"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "05_qc"), recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "geo_stage16/GEO_meta_ready_unique_gene_results_official_GSE300768.tsv.gz")
stopifnot(file.exists(input_file))
x <- fread(input_file, na.strings = c("NA", ""))

# Prespecified human evidence. Multiple contexts from the same GEO accession
# are combined first, so a study never contributes more than one effect per gene.
x[, context_id := paste(dataset, contrast, cell_line, sep = "|")]
x[, usable_effect := is.finite(log2FC) & is.finite(standard_error) & standard_error > 0 &
    !is.na(gene_symbol) & nzchar(gene_symbol)]

primary <- x[
  (dataset == "GSE200546" & contrast == "pH6.6_vs_pH7.4") |
  (dataset == "GSE116035" & contrast == "pH6.5_vs_pH7.4_adjusted_cell_line") |
  (dataset == "GSE152345" & contrast == "pH6.5_vs_pH7.6") |
  (dataset == "GSE29406" & contrast == "LA_vs_Ctrl") |
  (dataset == "GSE300768" & contrast == "pH6.4_vs_pH7.4") |
  (dataset == "GSE300765" & contrast %chin% c("acute_acidosis_vs_acute_neutral",
                                               "chronic_acidosis_vs_chronic_neutral")) |
  (dataset == "GSE9649" & contrast %chin% c("Acidosis_vs_Ctrl", "LA_vs_Ctrl"))
]
primary <- primary[usable_effect == TRUE]
secondary <- x[dataset == "GSE300768" & contrast == "pH6.4_vs_pH7.4" & usable_effect == TRUE]

context_audit <- unique(primary[, .(dataset, contrast, contrast_family, cell_line,
                                     analysis_method, platform, integration_role)])
context_audit[, meta_role := "primary"]
sec_audit <- unique(secondary[, .(dataset, contrast, contrast_family, cell_line,
                                   analysis_method, platform, integration_role)])
sec_audit[, meta_role := "secondary_sensitivity"]
context_audit <- rbind(context_audit, sec_audit, fill = TRUE)
setorder(context_audit, meta_role, dataset, contrast, cell_line)
fwrite(context_audit, file.path(out, "01_feasibility", "meta_context_inclusion_audit.tsv"), sep = "\t")

feasibility <- context_audit[, .(
  included_contexts = .N,
  included_datasets = uniqueN(dataset),
  RNAseq_contexts = sum(analysis_method == "DESeq2"),
  microarray_contexts = sum(analysis_method == "limma")
), by = meta_role]
fwrite(feasibility, file.path(out, "01_feasibility", "meta_feasibility_summary.tsv"), sep = "\t")

combine_within_dataset <- function(dt) {
  dt[, {
    w <- 1 / standard_error^2
    mu <- sum(w * log2FC) / sum(w)
    se <- sqrt(1 / sum(w))
    q <- sum(w * (log2FC - mu)^2)
    i2 <- if (.N > 1 && q > 0) max(0, (q - (.N - 1)) / q) * 100 else 0
    list(dataset_effect = mu, dataset_se = se, within_context_k = .N,
         within_Q = q, within_I2 = i2,
         within_direction_consistency = max(sum(log2FC > 0), sum(log2FC < 0)) / .N,
         contributing_contexts = paste(context_id, collapse = ";"))
  }, by = .(dataset, gene_symbol)]
}

primary_ds <- combine_within_dataset(primary)
secondary_ds <- combine_within_dataset(secondary)
fwrite(primary_ds, file.path(out, "02_dataset_level_effects", "primary_dataset_level_effects.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(secondary_ds, file.path(out, "02_dataset_level_effects", "secondary_GSE300768_effects.tsv.gz"),
       sep = "\t", compress = "gzip")

dl_core <- function(yi, sei) {
  ok <- is.finite(yi) & is.finite(sei) & sei > 0
  yi <- yi[ok]; sei <- sei[ok]
  k <- length(yi)
  if (!k) return(list(k = 0L))
  w <- 1 / sei^2
  fixed <- sum(w * yi) / sum(w)
  fixed_se <- sqrt(1 / sum(w))
  q <- sum(w * (yi - fixed)^2)
  df <- k - 1
  cval <- sum(w) - sum(w^2) / sum(w)
  tau2 <- if (k > 1 && cval > 0) max(0, (q - df) / cval) else 0
  wr <- 1 / (sei^2 + tau2)
  random <- sum(wr * yi) / sum(wr)
  random_se <- sqrt(1 / sum(wr))
  fixed_p <- 2 * pnorm(-abs(fixed / fixed_se))
  random_p <- 2 * pnorm(-abs(random / random_se))
  i2 <- if (k > 1 && q > 0) max(0, (q - df) / q) * 100 else 0
  dir_cons <- max(sum(yi > 0), sum(yi < 0)) / k
  list(k = k, fixed_effect = fixed, fixed_se = fixed_se, fixed_p = fixed_p,
       random_effect = random, random_se = random_se, random_p = random_p,
       Q = q, Q_p = if (k > 1) pchisq(q, df = df, lower.tail = FALSE) else NA_real_,
       I2 = i2, tau2 = tau2, direction_consistency = dir_cons,
       up_datasets = sum(yi > 0), down_datasets = sum(yi < 0))
}

dl_meta <- function(yi, sei) {
  ok <- is.finite(yi) & is.finite(sei) & sei > 0
  yi <- yi[ok]; sei <- sei[ok]
  ans <- dl_core(yi, sei)
  k <- ans$k
  loo_effect <- loo_p <- numeric(0)
  if (k >= 3) {
    for (j in seq_len(k)) {
      sub <- dl_core(yi[-j], sei[-j])
      loo_effect <- c(loo_effect, sub$random_effect)
      loo_p <- c(loo_p, sub$random_p)
    }
  }
  ans$loo_sign_stable <- if (k >= 3) all(sign(loo_effect) == sign(ans$random_effect)) else NA
  ans$loo_max_p <- if (k >= 3) max(loo_p) else NA_real_
  ans
}

run_meta <- function(ds, label) {
  ans <- ds[, as.list(dl_meta(dataset_effect, dataset_se)), by = gene_symbol]
  ans[, random_FDR := p.adjust(random_p, method = "BH")]
  ans[, fixed_FDR := p.adjust(fixed_p, method = "BH")]
  ans[, meta_direction := fifelse(random_effect > 0, "up_in_acidosis", "down_in_acidosis")]
  ans[, meta_supported := k >= 3]
  ans[, robust_candidate := k >= 3 & random_FDR < 0.05 & direction_consistency >= 0.70]
  ans[, strict_candidate := k >= 4 & random_FDR < 0.05 & direction_consistency >= 0.75 &
        loo_sign_stable == TRUE]
  setorder(ans, random_p)
  fwrite(ans, file.path(out, "03_human_meta_results", paste0(label, "_all_genes.tsv.gz")),
         sep = "\t", compress = "gzip")
  fwrite(ans[robust_candidate == TRUE],
         file.path(out, "03_human_meta_results", paste0(label, "_robust_FDR05.tsv")), sep = "\t")
  fwrite(ans[strict_candidate == TRUE],
         file.path(out, "03_human_meta_results", paste0(label, "_strict_FDR05.tsv")), sep = "\t")
  ans
}

# Broad acidosis is the primary analysis. Acute and chronic/lactic strata test
# whether signals survive a more homogeneous biological definition.
broad <- run_meta(primary_ds, "broad_acidosis_primary")

acute_contexts <- primary[contrast_family == "acute_acidosis"]
acute_ds <- combine_within_dataset(acute_contexts)
acute <- run_meta(acute_ds, "acute_acidosis_stratum")

chronic_contexts <- primary[contrast_family %chin% c("chronic_acidosis", "lactic_acidosis")]
chronic_ds <- combine_within_dataset(chronic_contexts)
chronic <- run_meta(chronic_ds, "chronic_lactic_acidosis_stratum")

broad_sensitivity_ds <- rbind(primary_ds, secondary_ds, fill = TRUE)
broad_secondary <- run_meta(broad_sensitivity_ds, "broad_acidosis_plus_GSE300768_sensitivity")

# Add stratum and secondary sensitivity evidence to the main table.
annot <- merge(broad, acute[, .(gene_symbol, acute_k = k, acute_effect = random_effect,
                                acute_p = random_p, acute_FDR = random_FDR)],
               by = "gene_symbol", all.x = TRUE)
annot <- merge(annot, chronic[, .(gene_symbol, chronic_lactic_k = k,
                                  chronic_lactic_effect = random_effect,
                                  chronic_lactic_p = random_p,
                                  chronic_lactic_FDR = random_FDR)],
               by = "gene_symbol", all.x = TRUE)
annot <- merge(annot, broad_secondary[, .(gene_symbol,
    secondary_inclusive_effect = random_effect,
    secondary_inclusive_p = random_p,
    secondary_inclusive_FDR = random_FDR)], by = "gene_symbol", all.x = TRUE)
annot[, secondary_direction_concordant :=
        sign(random_effect) == sign(secondary_inclusive_effect)]
setorder(annot, random_p)
fwrite(annot, file.path(out, "03_human_meta_results",
                        "broad_acidosis_primary_annotated_with_sensitivity.tsv.gz"),
       sep = "\t", compress = "gzip")

# Groessl CRISPR intersection using the versioned Ensembl release 113 audit.
crispr <- fread(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "s02/mouse_human_mapping_audit.tsv"), na.strings = c("", "NA"))
crispr[, support_n := fifelse(is.na(support_conditions) | support_conditions == "", 0L,
                              lengths(strsplit(support_conditions, ";", fixed = TRUE)))]
crispr[, human_symbol_key := toupper(human_symbol)]
annot[, human_symbol_key := toupper(gene_symbol)]

make_intersection <- function(meta_dt, allowed_status = NULL, label) {
  m <- copy(crispr)
  if (!is.null(allowed_status)) m <- m[mapping_status %chin% allowed_status]
  y <- merge(meta_dt, m, by = "human_symbol_key", allow.cartesian = TRUE)
  y <- unique(y, by = c("gene_symbol", "mouse_ensembl", "human_ensembl"))
  setorder(y, random_p, -support_n)
  fwrite(y, file.path(out, "04_crispr_intersections", paste0(label, ".tsv")), sep = "\t")
  y
}

robust_meta <- annot[robust_candidate == TRUE]
strict_meta <- annot[strict_candidate == TRUE]
inter_robust_1to1 <- make_intersection(robust_meta, "one_to_one",
                                       "robust_human_meta_x_CRISPR_one_to_one")
inter_strict_1to1 <- make_intersection(strict_meta, "one_to_one",
                                       "strict_human_meta_x_CRISPR_one_to_one")
inter_strict_all <- make_intersection(strict_meta, c("one_to_one", "one_to_many_or_mixed"),
                                      "strict_human_meta_x_CRISPR_all_mapped_sensitivity")
inter_robust_all <- make_intersection(robust_meta, c("one_to_one", "one_to_many_or_mixed"),
                                      "robust_human_meta_x_CRISPR_all_mapped_sensitivity")

intersection_summary <- data.table(
  set = c("robust_meta", "strict_meta", "robust_x_one_to_one", "strict_x_one_to_one",
          "robust_x_all_mapped", "strict_x_all_mapped"),
  n = c(nrow(robust_meta), nrow(strict_meta), uniqueN(inter_robust_1to1$gene_symbol),
        uniqueN(inter_strict_1to1$gene_symbol), uniqueN(inter_robust_all$gene_symbol),
        uniqueN(inter_strict_all$gene_symbol))
)
fwrite(intersection_summary, file.path(out, "04_crispr_intersections",
                                      "intersection_summary.tsv"), sep = "\t")

# Compact method/QC summaries.
meta_summary <- data.table(
  analysis = c("broad_primary", "acute_stratum", "chronic_lactic_stratum",
               "broad_plus_secondary"),
  datasets = c(uniqueN(primary_ds$dataset), uniqueN(acute_ds$dataset),
               uniqueN(chronic_ds$dataset), uniqueN(broad_sensitivity_ds$dataset)),
  tested_genes = c(nrow(broad), nrow(acute), nrow(chronic), nrow(broad_secondary)),
  genes_k3 = c(sum(broad$k >= 3), sum(acute$k >= 3), sum(chronic$k >= 3),
               sum(broad_secondary$k >= 3)),
  robust_FDR05 = c(sum(broad$robust_candidate), sum(acute$robust_candidate),
                   sum(chronic$robust_candidate), sum(broad_secondary$robust_candidate)),
  strict_FDR05 = c(sum(broad$strict_candidate), sum(acute$strict_candidate),
                   sum(chronic$strict_candidate), sum(broad_secondary$strict_candidate))
)
fwrite(meta_summary, file.path(out, "05_qc", "meta_analysis_summary.tsv"), sep = "\t")

pdf(file.path(out, "05_qc", "broad_meta_volcano.pdf"), width = 8, height = 7)
plot(broad$random_effect, -log10(pmax(broad$random_FDR, 1e-300)), pch = 16, cex = 0.35,
     col = ifelse(broad$strict_candidate, "#B2182B",
                  ifelse(broad$robust_candidate, "#2166AC", "#BDBDBD")),
     xlab = "Random-effects pooled log2 fold change",
     ylab = "-log10(BH FDR)", main = "Broad human acidosis meta-analysis")
abline(h = -log10(0.05), lty = 2)
dev.off()

writeLines(c(
  "Feasibility: PASS.",
  "Primary broad meta-analysis contains six independent GEO datasets.",
  "Multiple cell lines/contrasts inside one accession were inverse-variance combined before between-dataset meta-analysis.",
  "Between-dataset analysis reports fixed effects and DerSimonian-Laird random effects; random effects are primary.",
  "Robust candidate: k>=3 datasets, random-effects BH FDR<0.05, direction consistency>=0.70.",
  "Strict candidate: k>=4, random-effects BH FDR<0.05, direction consistency>=0.75, leave-one-dataset-out sign stable.",
  "GSE300768 is excluded from primary analysis and added only in the secondary-annotation sensitivity analysis.",
  "CRISPR primary intersection uses Ensembl release 113 one-to-one orthologs; all-mapped is sensitivity only.",
  "A statistical intersection is a candidate prioritization result, not functional validation."
), file.path(out, "METHOD_AND_RISK_CONTROL.txt"))

writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
cat("Completed. Broad robust:", nrow(robust_meta), "strict:", nrow(strict_meta),
    "strict x one-to-one:", uniqueN(inter_strict_1to1$gene_symbol), "\n")
