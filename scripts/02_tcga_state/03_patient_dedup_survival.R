args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript run_patient_dedup_and_survival.R <k2_survival.tsv> <uuid_mapping.tsv> <output_root>")
}

input_survival <- args[[1L]]
input_mapping <- args[[2L]]
output_root <- args[[3L]]

audit_dir <- file.path(output_root, "01_deduplication_audit")
patient_dir <- file.path(output_root, "02_patient_level_dataset")
survival_dir <- file.path(output_root, "03_survival_rerun")
session_dir <- file.path(output_root, "session_info")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(patient_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(survival_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(survival))

records <- read.delim(
  input_survival,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = ""
)
mapping <- read.delim(
  input_mapping,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = ""
)
names(mapping)[1L] <- sub("^\\ufeff", "", names(mapping)[1L])
if (!"file_uuid" %in% names(mapping)) names(mapping)[1L] <- "file_uuid"

stopifnot(nrow(records) == 647L)
stopifnot(!anyDuplicated(mapping[[1L]]))
mapping_index <- match(as.character(records[[1L]]), as.character(mapping[[1L]]))
records$patient_barcode_fixed <- as.character(mapping[[6L]][mapping_index])
if (anyNA(records$patient_barcode_fixed) || any(records$patient_barcode_fixed == "")) {
  stop("At least one expression UUID lacks a patient-barcode mapping")
}

records$sample_type_code <- substr(records$sample_barcode, 14L, 15L)
records$vial_code <- substr(records$sample_barcode, 16L, 16L)
records$sample_type_priority <- ifelse(records$sample_type_code == "01", 1L, 2L)
records$vial_priority <- match(records$vial_code, c(LETTERS, as.character(0:9)))
records$vial_priority[is.na(records$vial_priority)] <- 999L

ord <- order(
  records$patient_barcode_fixed,
  records$sample_type_priority,
  records$vial_priority,
  records$sample_barcode,
  records[[1L]]
)
records <- records[ord, , drop = FALSE]
records$selection_rank <- ave(
  seq_len(nrow(records)),
  records$patient_barcode_fixed,
  FUN = seq_along
)
records$retained_patient_record <- records$selection_rank == 1L
records$selection_reason <- ifelse(
  records$retained_patient_record,
  "retained_by_fixed_barcode_then_uuid_rule",
  "excluded_additional_primary_tumor_record"
)

patient_counts <- table(records$patient_barcode_fixed)
records$records_for_patient <- unname(patient_counts[records$patient_barcode_fixed])
cluster_counts <- tapply(records$cluster, records$patient_barcode_fixed, function(x) length(unique(x)))
records$discordant_cluster_within_patient <- unname(cluster_counts[records$patient_barcode_fixed]) > 1L

patients <- records[records$retained_patient_record, , drop = FALSE]
stopifnot(nrow(patients) == 624L)
stopifnot(!anyDuplicated(patients$patient_barcode_fixed))

write.table(
  records,
  file.path(audit_dir, "TCGA_647_record_fixed_selection_audit.tsv"),
  sep = "\t", row.names = FALSE, quote = TRUE, na = ""
)
write.table(
  records[records$records_for_patient > 1L, , drop = FALSE],
  file.path(audit_dir, "TCGA_13_duplicate_patient_record_audit.tsv"),
  sep = "\t", row.names = FALSE, quote = TRUE, na = ""
)
write.table(
  records[!records$retained_patient_record, , drop = FALSE],
  file.path(audit_dir, "TCGA_23_excluded_records.tsv"),
  sep = "\t", row.names = FALSE, quote = TRUE, na = ""
)
write.table(
  patients,
  file.path(patient_dir, "TCGA_COADREAD_K2_patient_level_624.tsv"),
  sep = "\t", row.names = FALSE, quote = TRUE, na = ""
)

rule_text <- c(
  "Fixed patient-level record selection rule",
  "1. Group Primary Tumor expression records by the mapped 12-character TCGA patient barcode.",
  "2. Prefer sample-type code 01.",
  "3. Prefer vial codes alphabetically (A before B before C, followed by other codes).",
  "4. Break remaining ties by sample barcode and then file UUID in ascending lexical order.",
  "5. Do not use C1/C2 label, survival outcome, CMS, MSI, effect size, or statistical significance.",
  "6. Apply endpoint-specific complete-case filtering only after one record per patient has been retained."
)
writeLines(rule_text, file.path(audit_dir, "fixed_deduplication_rule.txt"), useBytes = TRUE)

duplicate_ids <- names(patient_counts)[patient_counts > 1L]
summary_rows <- data.frame(
  metric = c(
    "starting_primary_tumor_records",
    "unique_sample_barcodes",
    "unique_patients",
    "patients_with_multiple_records",
    "excluded_additional_records",
    "duplicate_patients_with_discordant_K2_labels"
  ),
  value = c(
    nrow(records),
    length(unique(records$sample_barcode)),
    length(unique(records$patient_barcode_fixed)),
    length(duplicate_ids),
    sum(!records$retained_patient_record),
    sum(cluster_counts[duplicate_ids] > 1L)
  )
)
write.table(
  summary_rows,
  file.path(audit_dir, "deduplication_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = TRUE
)

endpoints <- c("OS", "PFI", "DFI", "DSS")
results <- vector("list", length(endpoints))
names(results) <- endpoints
km_fits <- vector("list", length(endpoints))
names(km_fits) <- endpoints

for (endpoint in endpoints) {
  time_col <- paste0(endpoint, ".time")
  analysis <- patients[
    !is.na(patients[[endpoint]]) & patients[[endpoint]] != "" &
      !is.na(patients[[time_col]]) & patients[[time_col]] != "",
    , drop = FALSE
  ]
  analysis$event <- as.integer(analysis[[endpoint]])
  analysis$time_days <- as.numeric(analysis[[time_col]])
  analysis$cluster_factor <- factor(analysis$cluster, levels = c(1, 2), labels = c("C1", "C2"))
  analysis <- analysis[!is.na(analysis$event) & !is.na(analysis$time_days) & !is.na(analysis$cluster_factor), , drop = FALSE]

  fit <- survfit(Surv(time_days, event) ~ cluster_factor, data = analysis)
  test <- survdiff(Surv(time_days, event) ~ cluster_factor, data = analysis)
  logrank_p <- pchisq(test$chisq, df = length(test$n) - 1L, lower.tail = FALSE)
  cox <- coxph(Surv(time_days, event) ~ cluster_factor, data = analysis)
  cox_summary <- summary(cox)

  results[[endpoint]] <- data.frame(
    endpoint = endpoint,
    n_patients = nrow(analysis),
    n_C1 = sum(analysis$cluster_factor == "C1"),
    n_C2 = sum(analysis$cluster_factor == "C2"),
    events_total = sum(analysis$event == 1L),
    events_C1 = sum(analysis$event == 1L & analysis$cluster_factor == "C1"),
    events_C2 = sum(analysis$event == 1L & analysis$cluster_factor == "C2"),
    logrank_chisq = unname(test$chisq),
    logrank_p = logrank_p,
    cox_HR_C2_vs_C1 = unname(cox_summary$coefficients[1L, "exp(coef)"]),
    cox_CI_low = unname(cox_summary$conf.int[1L, "lower .95"]),
    cox_CI_high = unname(cox_summary$conf.int[1L, "upper .95"]),
    cox_p = unname(cox_summary$coefficients[1L, "Pr(>|z|)"])
  )
  km_fits[[endpoint]] <- fit
  write.table(
    analysis,
    file.path(survival_dir, paste0(endpoint, "_patient_level_analysis_dataset.tsv")),
    sep = "\t", row.names = FALSE, quote = TRUE, na = ""
  )
}

results <- do.call(rbind, results)
write.table(
  results,
  file.path(survival_dir, "K2_four_endpoint_patient_level_survival.tsv"),
  sep = "\t", row.names = FALSE, quote = TRUE
)

old_results <- data.frame(
  endpoint = c("OS", "PFI", "DFI", "DSS"),
  old_record_n = c(643L, 643L, 247L, 621L),
  old_record_events = c(129L, 164L, 29L, 78L),
  old_record_logrank_p = c(0.197502203920908, 0.383124696280239, 0.950128962870498, 0.945548465307597)
)
comparison <- merge(old_results, results, by = "endpoint", sort = FALSE)
comparison <- comparison[match(endpoints, comparison$endpoint), , drop = FALSE]
write.table(
  comparison,
  file.path(survival_dir, "record_level_vs_patient_level_survival_comparison.tsv"),
  sep = "\t", row.names = FALSE, quote = TRUE
)

png(file.path(survival_dir, "K2_four_endpoint_patient_level_KM.png"), width = 2400, height = 2000, res = 300)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3.0, 1.0))
for (endpoint in endpoints) {
  row <- results[results$endpoint == endpoint, ]
  plot(
    km_fits[[endpoint]],
    col = c("#2C7FB8", "#D95F59"), lwd = 2,
    xlab = "Time (days)", ylab = "Survival probability",
    main = sprintf("%s (n=%d, log-rank P=%.3g)", endpoint, row$n_patients, row$logrank_p),
    mark.time = TRUE
  )
  legend("bottomleft", legend = c("C1", "C2"), col = c("#2C7FB8", "#D95F59"), lwd = 2, bty = "n")
}
dev.off()

writeLines(capture.output(sessionInfo()), file.path(session_dir, "R_sessionInfo_patient_dedup_survival.txt"))

input_manifest <- data.frame(
  role = c("K2 survival input", "UUID-barcode mapping"),
  path = c(normalizePath(input_survival, winslash = "/"), normalizePath(input_mapping, winslash = "/")),
  md5 = unname(tools::md5sum(c(input_survival, input_mapping)))
)
write.table(
  input_manifest,
  file.path(session_dir, "input_file_manifest.tsv"),
  sep = "\t", row.names = FALSE, quote = TRUE
)

message("Completed: 647 records -> ", nrow(patients), " patients")
print(results)
