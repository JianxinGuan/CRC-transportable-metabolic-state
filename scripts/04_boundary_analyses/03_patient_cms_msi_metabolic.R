args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop("Usage: Rscript script patient.tsv acidosis.gz msi.tsv cms.tsv output_root session_dir")
patient_file <- args[[1L]]; acidosis_file <- args[[2L]]; msi_file <- args[[3L]]; cms_file <- args[[4L]]; out_root <- args[[5L]]; session_dir <- args[[6L]]
out <- file.path(out_root, "04_CMS_MSI_rerun")
dir.create(out, recursive = TRUE, showWarnings = FALSE); dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(stats))

read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
patients <- read_tsv(patient_file)
marker <- read.delim(gzfile(acidosis_file), check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
msi <- read_tsv(msi_file)
cms <- read_tsv(cms_file)
names(msi)[1L] <- sub("^\\ufeff", "", names(msi)[1L])

if (!"sample_id" %in% names(marker)) stop("marker table lacks sample_id")
if (!"sample_id" %in% names(patients)) stop("patient table lacks sample_id")
if (!"patient_barcode" %in% names(patients)) {
  if ("patient_barcode_fixed" %in% names(patients)) patients$patient_barcode <- patients$patient_barcode_fixed else stop("patient table lacks patient_barcode")
}

## CMSclassifier sample-level result: first column is UUID and predicted CMS is last column.
names(cms)[1L] <- "sample_id"
cms_col <- if ("CMS" %in% names(cms)) "CMS" else if ("RF.details.RF.predictedCMS" %in% names(cms)) "RF.details.RF.predictedCMS" else tail(names(cms), 1L)
cms_small <- cms[, c("sample_id", cms_col), drop = FALSE]
names(cms_small)[2L] <- "CMS"

marker_keep <- marker[match(patients$sample_id, marker$sample_id), , drop = FALSE]
msi_by_patient <- msi[order(msi$patient_barcode, msi$sample_barcode, msi[[1L]]), , drop = FALSE]
msi_by_patient <- msi_by_patient[!duplicated(msi_by_patient$patient_barcode), , drop = FALSE]
msi_keep <- msi_by_patient[match(patients$patient_barcode, msi_by_patient$patient_barcode), , drop = FALSE]
cms_keep <- cms_small[match(patients$sample_id, cms_small$sample_id), , drop = FALSE]

dat <- patients
dat$CMS <- ifelse(!is.na(marker_keep$CMS) & marker_keep$CMS != "", as.character(marker_keep$CMS), as.character(cms_keep$CMS))
dat$CMS <- as.character(dat$CMS)
for (nm in c("carbonic_anhydrase", "acid_extrusion_transport", "lactate_axis", "ca9_ca12", "acidosis_composite", "CA9", "CA12", "acidosis_signed_score", "hypoxia_score", "MicroenvironmentScore", "StromaScore", "ImmuneScore", "metabolic_intensity", "proliferation", "microenv_proxy")) {
  if (nm %in% names(marker_keep)) dat[[nm]] <- suppressWarnings(as.numeric(marker_keep[[nm]]))
}
dat$MANTIS <- suppressWarnings(as.numeric(msi_keep$MSI_SCORE_MANTIS))
dat$MSIsensor <- suppressWarnings(as.numeric(msi_keep$MSI_SENSOR_SCORE))
dat$MANTIS_MSIH_like <- ifelse(!is.na(dat$MANTIS), ifelse(dat$MANTIS >= 0.4, "MSI-H-like", "MSS-like"), NA)
dat$MSIsensor_MSIH_like <- ifelse(!is.na(dat$MSIsensor), ifelse(dat$MSIsensor >= 0.5, "MSI-H-like", "MSS-like"), NA)
dat$state <- factor(dat$cluster, levels = c(1, 2), labels = c("C1", "C2"))
dat$CMS <- factor(dat$CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
write.table(dat, file.path(out, "TCGA_624_patient_level_CMS_MSI_metabolic_merged.tsv"), sep = "\t", row.names = FALSE, quote = TRUE, na = "")

## K2-CMS global association and effect size.
cms_dat <- dat[!is.na(dat$CMS), , drop = FALSE]
tab_cms <- table(cms_dat$state, cms_dat$CMS)
chis <- suppressWarnings(chisq.test(tab_cms, correct = FALSE))
fisher <- fisher.test(tab_cms, workspace = 2e8)
cramer_v <- sqrt(as.numeric(chis$statistic) / (sum(tab_cms) * min(nrow(tab_cms)-1, ncol(tab_cms)-1)))
write.table(as.data.frame.matrix(tab_cms), file.path(out, "K2_CMS_contingency_624_patient_level.tsv"), sep = "\t", quote = TRUE)
write.table(data.frame(n = sum(tab_cms), chi_square = unname(chis$statistic), chi_df = unname(chis$parameter), chi_p = chis$p.value, fisher_p = fisher$p.value, cramers_V = cramer_v, min_expected_count = min(chis$expected)), file.path(out, "K2_CMS_global_test_624_patient_level.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)

## Within-CMS Wilcoxon comparisons for acidosis composite and metabolic intensity.
themes <- c("acidosis_composite", "metabolic_intensity")
within <- list(); idx <- 1L
for (theme in themes) for (cc in levels(dat$CMS)) {
  x <- dat[!is.na(dat$CMS) & dat$CMS == cc & !is.na(dat[[theme]]) & !is.na(dat$state), , drop = FALSE]
  if (nrow(x) < 2L || length(unique(x$state)) < 2L) next
  wt <- wilcox.test(x[[theme]] ~ x$state, exact = FALSE)
  means <- tapply(x[[theme]], x$state, mean); meds <- tapply(x[[theme]], x$state, median)
  within[[idx]] <- data.frame(outcome = theme, CMS = cc, n = nrow(x), C1_n = sum(x$state == "C1"), C2_n = sum(x$state == "C2"), mean_C1 = unname(means["C1"]), mean_C2 = unname(means["C2"]), median_C1 = unname(meds["C1"]), median_C2 = unname(meds["C2"]), delta_C2_minus_C1 = unname(medians <- meds["C2"] - meds["C1"]), wilcox_p = wt$p.value)
  idx <- idx + 1L
}
within <- do.call(rbind, within); within$BH_FDR <- p.adjust(within$wilcox_p, method = "BH")
write.table(within, file.path(out, "within_CMS_metabolic_effects_624_patient_level.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)

## MSI continuous Wilcoxon and binary Fisher tests.
msi_rows <- list(); idx <- 1L
for (nm in c("MANTIS", "MSIsensor")) {
  x <- dat[!is.na(dat[[nm]]) & !is.na(dat$state), , drop = FALSE]
  wt <- wilcox.test(x[[nm]] ~ x$state, exact = FALSE)
  group_med <- tapply(x[[nm]], x$state, median)
  msi_rows[[idx]] <- data.frame(metric = nm, analysis = "continuous_Wilcoxon", n = nrow(x), C1_n = sum(x$state == "C1"), C2_n = sum(x$state == "C2"), C1_median = unname(group_med["C1"]), C2_median = unname(group_med["C2"]), effect_C2_minus_C1 = unname(group_med["C2"] - group_med["C1"]), p = wt$p.value, FDR = NA_real_); idx <- idx + 1L
  catvar <- paste0(nm, "_MSIH_like"); y <- dat[!is.na(dat[[catvar]]) & !is.na(dat$state), , drop = FALSE]; ft <- fisher.test(table(y$state, y[[catvar]])); ytab <- as.data.frame.matrix(table(y$state, y[[catvar]])); write.table(ytab, file.path(out, paste0(nm, "_MSIH_like_contingency_624.tsv")), sep = "\t", quote = TRUE)
  msi_rows[[idx]] <- data.frame(metric = nm, analysis = "binary_Fisher", n = nrow(y), C1_n = sum(y$state == "C1"), C2_n = sum(y$state == "C2"), C1_median = NA, C2_median = NA, effect_C2_minus_C1 = NA, p = ft$p.value, FDR = NA_real_); idx <- idx + 1L
}
msi_out <- do.call(rbind, msi_rows); msi_out$FDR <- p.adjust(msi_out$p, method = "BH")
write.table(msi_out, file.path(out, "MSI_continuous_and_binary_K2_tests_624_patient_level.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
writeLines(c("MANTIS MSI-H-like threshold: >= 0.4", "MSIsensor MSI-H-like threshold: >= 0.5", "Thresholds retained from the prior frozen MSI audit; verify against the original PanCanAtlas documentation before publication."), file.path(out, "MSI_thresholds_and_limitations.txt"))

## CMS incremental linear models using available prespecified covariates.
model_rows <- list(); idx <- 1L
for (theme in c("acidosis_composite", "metabolic_intensity")) {
  covars <- c("proliferation", "hypoxia_score", "microenv_proxy", "CMS", "primary_site", "stage_group")
  use <- c(theme, "state", covars); use <- use[use %in% names(dat)]
  x <- dat[complete.cases(dat[, use, drop = FALSE]), , drop = FALSE]
  if (nrow(x) < 30L) next
  f0 <- as.formula(paste(theme, "~", paste(setdiff(use, c(theme, "state")), collapse = " + ")))
  f1 <- as.formula(paste(theme, "~ state +", paste(setdiff(use, c(theme, "state")), collapse = " + ")))
  m0 <- lm(f0, data = x); m1 <- lm(f1, data = x); sm <- summary(m1); an <- anova(m0, m1)
  coefrow <- sm$coefficients[grep("^state", rownames(sm$coefficients))[1L], , drop = FALSE]
  model_rows[[idx]] <- data.frame(outcome = theme, analysis = "CMS_incremental_available_covariates", n = nrow(x), state_beta = coefrow[1,"Estimate"], SE = coefrow[1,"Std. Error"], CI_low = confint(m1)[grep("^state", rownames(confint(m1)))[1L],1], CI_high = confint(m1)[grep("^state", rownames(confint(m1)))[1L],2], state_p = coefrow[1,"Pr(>|t|)"], LRT_F = an$F[2], LRT_p = an$`Pr(>F)`[2], adjusted_R2_base = summary(m0)$adj.r.squared, adjusted_R2_full = sm$adj.r.squared, adjusted_R2_increment = sm$adj.r.squared - summary(m0)$adj.r.squared)
  idx <- idx + 1L
}
model_out <- do.call(rbind, model_rows); model_out$FDR <- p.adjust(model_out$state_p, method = "BH")
write.table(model_out, file.path(out, "CMS_incremental_models_624_patient_level.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)

writeLines(capture.output(sessionInfo()), file.path(session_dir, "R_sessionInfo_CMS_MSI_metabolic_624.txt"))
write.table(data.frame(file = c(patient_file, acidosis_file, msi_file, cms_file), md5 = unname(tools::md5sum(c(patient_file, acidosis_file, msi_file, cms_file)))), file.path(session_dir, "CMS_MSI_input_manifest.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
message("Completed patient-level CMS/MSI/metabolic rerun: ", nrow(dat), " patients; CMS evaluable=", nrow(cms_dat))
