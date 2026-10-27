invisible(Sys.setlocale("LC_CTYPE", "Chinese"))
out <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "geo_stage15")
full_file <- file.path(out, "unified_results", "GEO_all_contrasts_unified.tsv.gz")
x <- read.delim(gzfile(full_file), check.names = FALSE, stringsAsFactors = FALSE)

valid <- x$gene_symbol_valid
context <- paste(x$dataset, x$contrast, x$cell_line, sep = "|")
key <- paste(context, x$gene_symbol, sep = "|")
direction_mismatch <- sum(
  (x$log2FC > 0 & x$direction != "up_in_treatment") |
    (x$log2FC < 0 & x$direction != "down_in_treatment"), na.rm = TRUE)
duplicate_valid_keys <- sum(duplicated(key) & valid)

# A Meta-ready table requires one symbol per analysis context. When multiple
# source IDs map to one symbol, retain the row with smallest P value, then
# largest absolute effect. The complete non-deduplicated table remains intact.
z <- x[valid & !is.na(x$P_Value) & !is.na(x$log2FC), ]
z$context_key <- paste(z$dataset, z$contrast, z$cell_line, sep = "|")
z$gene_key <- paste(z$context_key, z$gene_symbol, sep = "|")
ord <- order(z$gene_key, z$P_Value, -abs(z$log2FC), z$gene_id)
z <- z[ord, ]
z$mapping_duplicate_count <- ave(z$gene_key, z$gene_key, FUN = length)
meta_ready <- z[!duplicated(z$gene_key), ]
meta_ready$gene_key <- NULL

con <- gzfile(file.path(out, "unified_results", "GEO_meta_ready_unique_gene_results.tsv.gz"), "wt")
write.table(meta_ready, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
close(con)

sig_check <- read.delim(file.path(out, "unified_results", "GEO_nominal_logFC2_P005_unified.tsv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
expected_sig <- sum(x$pass_abs_log2FC2_P005 & valid)

missing_by_dataset <- aggregate(!valid, list(dataset = x$dataset), sum)
colnames(missing_by_dataset)[2] <- "rows_missing_gene_symbol"
write.table(missing_by_dataset,
            file.path(out, "unified_results", "gene_symbol_missingness_by_dataset.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

checks <- data.frame(
  check = c("full_rows", "valid_gene_symbol_rows", "missing_gene_symbol_rows",
            "analysis_contexts", "duplicate_valid_context_gene_keys",
            "meta_ready_rows", "direction_mismatches", "expected_significant_rows",
            "written_significant_rows", "significant_count_matches",
            "audit_dataset_rows"),
  value = as.character(c(
    nrow(x), sum(valid), sum(!valid), length(unique(context)), duplicate_valid_keys,
    nrow(meta_ready), direction_mismatch, expected_sig, nrow(sig_check),
    expected_sig == nrow(sig_check),
    nrow(read.delim(file.path(out, "final_cohort_eligibility_audit.tsv")))
  )), stringsAsFactors = FALSE
)
write.table(checks, file.path(out, "VALIDATION_REPORT.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

stopifnot(direction_mismatch == 0L)
stopifnot(expected_sig == nrow(sig_check))
stopifnot(length(unique(context)) == 20L)
stopifnot(!anyDuplicated(paste(meta_ready$context_key, meta_ready$gene_symbol, sep = "|")))
stopifnot(nrow(read.delim(file.path(out, "final_cohort_eligibility_audit.tsv"))) == 9L)

cat("Validation passed. Full rows:", nrow(x), "Meta-ready rows:", nrow(meta_ready), "\n")
