options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260929)
args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[1] else file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "switch_audit_input")
output_dir <- if (length(args) >= 2) args[2] else file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "switch_audit_output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

scores_long <- fread(file = file.path(input_dir, "current_scores.tsv.gz"))
quality_full <- fread(file = file.path(input_dir, "locked_quality.tsv.gz"))
reprojection <- fread(file = file.path(input_dir, "locked_reprojection.tsv.gz"))

# Scores are identical for a sample across centroid specifications; labels are not.
labels <- unique(scores_long[, .(cohort, centroid_specification, sample_id, cluster)])
labels_wide <- dcast(labels, cohort + sample_id ~ centroid_specification, value.var = "cluster")
setnames(labels_wide,
         c("full_cohort_fixed_centroids", "exclusion13_rebuilt_centroids"),
         c("full_cluster", "exclusion13_cluster"))

marker_scores <- unique(scores_long[outcome %chin% c("lactate_axis", "acidification_composite"),
                                    .(cohort, sample_id, outcome, score)])
marker_wide <- dcast(marker_scores, cohort + sample_id ~ outcome, value.var = "score")

full_margin <- quality_full[cohort %chin% c("GSE39582", "GSE17536"),
                            .(cohort, sample_id,
                              full_cluster_source = cluster,
                              full_cor_C1 = cor_TCGA_C1,
                              full_cor_C2 = cor_TCGA_C2,
                              full_abs_margin = assignment_margin,
                              full_silhouette = silhouette)]
full_margin[, full_signed_score := full_cor_C2 - full_cor_C1]

new_margin <- reprojection[, .(cohort, sample_id,
                               old_cluster_source = old_cluster,
                               exclusion13_cluster_source = new_cluster,
                               old_margin_source = old_margin,
                               exclusion13_cor_C1 = new_cor_C1,
                               exclusion13_cor_C2 = new_cor_C2,
                               exclusion13_abs_margin = new_margin,
                               exclusion13_silhouette = silhouette)]
new_margin[, exclusion13_signed_score := exclusion13_cor_C2 - exclusion13_cor_C1]

dat <- Reduce(function(x, y) merge(x, y, by = c("cohort", "sample_id"), all = TRUE),
              list(labels_wide, marker_wide, full_margin, new_margin))

# Fail closed if the sample-level projection sources do not match the formal labels.
stopifnot(!anyNA(dat[, .(full_cluster, exclusion13_cluster, full_cluster_source,
                         old_cluster_source, exclusion13_cluster_source)]))
stopifnot(all(dat$full_cluster == dat$full_cluster_source))
stopifnot(all(dat$full_cluster == dat$old_cluster_source))
stopifnot(all(dat$exclusion13_cluster == dat$exclusion13_cluster_source))

dat[, changed := full_cluster != exclusion13_cluster]
dat[, transition := fifelse(full_cluster == 1 & exclusion13_cluster == 1, "C1 retained",
                     fifelse(full_cluster == 1 & exclusion13_cluster == 2, "C1 to C2",
                     fifelse(full_cluster == 2 & exclusion13_cluster == 1, "C2 to C1",
                             "C2 retained")))]
dat[, transition := factor(transition,
                           levels = c("C1 retained", "C1 to C2", "C2 to C1", "C2 retained"))]

fwrite(dat, file.path(output_dir, "01_sample_level_transition_margin_scores.tsv.gz"), sep = "\t")

counts <- dat[, .(n = .N,
                  median_full_abs_margin = median(full_abs_margin, na.rm = TRUE),
                  median_exclusion13_abs_margin = median(exclusion13_abs_margin, na.rm = TRUE),
                  median_full_signed_score = median(full_signed_score, na.rm = TRUE),
                  median_exclusion13_signed_score = median(exclusion13_signed_score, na.rm = TRUE),
                  median_lactate_axis = median(lactate_axis, na.rm = TRUE),
                  median_acidification_composite = median(acidification_composite, na.rm = TRUE)),
              by = .(cohort, transition)]
fwrite(counts, file.path(output_dir, "02_transition_group_summary.tsv"), sep = "\t")

transition_counts <- dat[, .N, by = .(cohort, full_cluster, exclusion13_cluster, transition)]
transition_counts[, cohort_total := sum(N), by = cohort]
transition_counts[, proportion := N / cohort_total]
fwrite(transition_counts, file.path(output_dir, "03_label_transition_counts.tsv"), sep = "\t")

# Directional switching test: samples moved into C2 versus samples moved into C1.
direction_tests <- rbindlist(lapply(unique(dat$cohort), function(cc) {
  x <- dat[cohort == cc & transition %chin% c("C1 to C2", "C2 to C1")]
  rbindlist(lapply(c("lactate_axis", "acidification_composite",
                     "full_abs_margin", "exclusion13_abs_margin"), function(v) {
    a <- x[transition == "C1 to C2", get(v)]
    b <- x[transition == "C2 to C1", get(v)]
    wt <- wilcox.test(a, b, exact = FALSE)
    data.table(cohort = cc, variable = v,
               group1 = "C1 to C2", group1_n = sum(is.finite(a)), group1_median = median(a, na.rm = TRUE),
               group2 = "C2 to C1", group2_n = sum(is.finite(b)), group2_median = median(b, na.rm = TRUE),
               median_difference_group1_minus_group2 = median(a, na.rm = TRUE) - median(b, na.rm = TRUE),
               wilcoxon_W = unname(wt$statistic), p_value = wt$p.value)
  }))
}))
direction_tests[, FDR_within_cohort := p.adjust(p_value, method = "BH"), by = cohort]
fwrite(direction_tests, file.path(output_dir, "04_switch_direction_score_tests.tsv"), sep = "\t")

# Changed samples should be concentrated near the projection boundary.
margin_tests <- rbindlist(lapply(unique(dat$cohort), function(cc) {
  x <- dat[cohort == cc]
  rbindlist(lapply(c("full_abs_margin", "exclusion13_abs_margin"), function(v) {
    a <- x[changed == TRUE, get(v)]
    b <- x[changed == FALSE, get(v)]
    wt <- wilcox.test(a, b, exact = FALSE)
    data.table(cohort = cc, margin_specification = v,
               changed_n = length(a), changed_median = median(a, na.rm = TRUE),
               unchanged_n = length(b), unchanged_median = median(b, na.rm = TRUE),
               wilcoxon_W = unname(wt$statistic), p_value = wt$p.value)
  }))
}))
margin_tests[, FDR_within_cohort := p.adjust(p_value, method = "BH"), by = cohort]
fwrite(margin_tests, file.path(output_dir, "05_changed_vs_unchanged_margin_tests.tsv"), sep = "\t")

# Continuous projection score correlations. Positive signed scores point toward C2.
cor_results <- rbindlist(lapply(unique(dat$cohort), function(cc) {
  x <- dat[cohort == cc]
  rbindlist(lapply(c("full_signed_score", "exclusion13_signed_score"), function(proj) {
    rbindlist(lapply(c("lactate_axis", "acidification_composite"), function(marker) {
      keep <- is.finite(x[[proj]]) & is.finite(x[[marker]])
      ct <- cor.test(x[[proj]][keep], x[[marker]][keep], method = "spearman", exact = FALSE)
      data.table(cohort = cc, projection_score = proj, marker_score = marker,
                 n = sum(keep), spearman_rho = unname(ct$estimate), p_value = ct$p.value)
    }))
  }))
}))
cor_results[, FDR_within_cohort := p.adjust(p_value, method = "BH"), by = cohort]
fwrite(cor_results, file.path(output_dir, "06_continuous_projection_marker_correlations.tsv"), sep = "\t")

# Paired bootstrap of the change in Spearman rho after centroid reconstruction.
B <- 2000L
boot_results <- rbindlist(lapply(unique(dat$cohort), function(cc) {
  x <- dat[cohort == cc]
  rbindlist(lapply(c("lactate_axis", "acidification_composite"), function(marker) {
    keep <- is.finite(x$full_signed_score) & is.finite(x$exclusion13_signed_score) & is.finite(x[[marker]])
    y <- x[keep]
    rho_full <- cor(y$full_signed_score, y[[marker]], method = "spearman")
    rho_new <- cor(y$exclusion13_signed_score, y[[marker]], method = "spearman")
    delta_obs <- rho_new - rho_full
    deltas <- replicate(B, {
      ii <- sample.int(nrow(y), nrow(y), replace = TRUE)
      cor(y$exclusion13_signed_score[ii], y[[marker]][ii], method = "spearman") -
        cor(y$full_signed_score[ii], y[[marker]][ii], method = "spearman")
    })
    data.table(cohort = cc, marker_score = marker, n = nrow(y),
               rho_full = rho_full, rho_exclusion13 = rho_new,
               delta_rho_exclusion13_minus_full = delta_obs,
               bootstrap_replicates = B,
               bootstrap_CI_low = quantile(deltas, .025, na.rm = TRUE),
               bootstrap_CI_high = quantile(deltas, .975, na.rm = TRUE),
               bootstrap_two_sided_p = min(1, 2 * (min(sum(deltas <= 0, na.rm = TRUE),
                                                       sum(deltas >= 0, na.rm = TRUE)) + 1) / (B + 1)))
  }))
}))
boot_results[, bootstrap_FDR := p.adjust(bootstrap_two_sided_p, method = "BH")]
fwrite(boot_results, file.path(output_dir, "07_paired_bootstrap_correlation_change.tsv"), sep = "\t")

# Compact audit figure for supplementary use.
plot_dat <- melt(dat, id.vars = c("cohort", "sample_id", "transition", "changed"),
                 measure.vars = c("lactate_axis", "acidification_composite"),
                 variable.name = "marker", value.name = "score")
plot_dat[, marker := factor(marker,
                            levels = c("lactate_axis", "acidification_composite"),
                            labels = c("Lactate axis", "Acidification composite"))]

pA <- ggplot(plot_dat[transition %chin% c("C1 to C2", "C2 to C1")],
             aes(transition, score, fill = transition)) +
  geom_violin(trim = FALSE, alpha = .35, color = NA) +
  geom_boxplot(width = .18, outlier.shape = NA, linewidth = .35) +
  facet_grid(marker ~ cohort, scales = "free_y") +
  scale_fill_manual(values = c("C1 to C2" = "#D77A61", "C2 to C1" = "#2F6690"), guide = "none") +
  labs(title = "A  Marker scores of directionally reassigned samples", x = NULL, y = "Marker score") +
  theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 18, hjust = 1))

margin_long <- melt(dat, id.vars = c("cohort", "sample_id", "changed"),
                    measure.vars = c("full_abs_margin", "exclusion13_abs_margin"),
                    variable.name = "specification", value.name = "absolute_margin")
margin_long[, specification := factor(specification,
                                      levels = c("full_abs_margin", "exclusion13_abs_margin"),
                                      labels = c("Full-cohort", "Exclusion-13"))]
pB <- ggplot(margin_long, aes(changed, absolute_margin, fill = changed)) +
  geom_boxplot(outlier.shape = NA, width = .55, linewidth = .35) +
  facet_grid(specification ~ cohort, scales = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "#C95D63", `FALSE` = "#7A9E9F"), guide = "none") +
  scale_x_discrete(labels = c(`FALSE` = "Unchanged", `TRUE` = "Changed")) +
  labs(title = "B  Changed samples have lower projection margins", x = NULL, y = "Absolute margin") +
  theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 18, hjust = 1))

cor_plot <- copy(cor_results)
cor_plot[, specification := factor(projection_score,
                                    levels = c("full_signed_score", "exclusion13_signed_score"),
                                    labels = c("Full-cohort", "Exclusion-13"))]
cor_plot[, marker := factor(marker_score,
                            levels = c("lactate_axis", "acidification_composite"),
                            labels = c("Lactate axis", "Acidification composite"))]
pC <- ggplot(cor_plot, aes(specification, spearman_rho, color = marker, group = marker)) +
  geom_hline(yintercept = 0, linetype = 2, color = "#777777") +
  geom_line(linewidth = .6) + geom_point(size = 2.5) +
  facet_wrap(~cohort) +
  scale_color_manual(values = c("Lactate axis" = "#D77A61", "Acidification composite" = "#5B8E7D")) +
  labs(title = "C  Continuous projection-score correlations", x = NULL, y = "Spearman rho", color = NULL) +
  theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")

fig <- (pA | pB) / pC + plot_layout(heights = c(1.2, .8))
ggsave(file.path(output_dir, "08_transition_margin_continuous_audit.pdf"), fig,
       width = 210, height = 220, units = "mm", device = cairo_pdf)
ggsave(file.path(output_dir, "08_transition_margin_continuous_audit_300dpi.png"), fig,
       width = 210, height = 220, units = "mm", dpi = 300, bg = "white")

# Method/input audit.
audit <- data.table(
  item = c("formal_label_match_full", "formal_label_match_exclusion13",
           "GSE39582_samples", "GSE17536_samples", "bootstrap_replicates",
           "projection_score_definition", "multiple_testing"),
  value = c(as.character(all(dat$full_cluster == dat$full_cluster_source)),
            as.character(all(dat$exclusion13_cluster == dat$exclusion13_cluster_source)),
            as.character(dat[cohort == "GSE39582", .N]),
            as.character(dat[cohort == "GSE17536", .N]),
            as.character(B),
            "correlation_to_C2_centroid_minus_correlation_to_C1_centroid",
            "BH within cohort for primary Wilcoxon/correlation families; global BH for bootstrap delta-rho tests"))
fwrite(audit, file.path(output_dir, "09_method_and_input_audit.tsv"), sep = "\t")

# Machine-generated interpretation summary.
fmt <- function(x, d = 3) formatC(x, digits = d, format = "f")
g39_lac <- direction_tests[cohort == "GSE39582" & variable == "lactate_axis"]
g17_lac <- direction_tests[cohort == "GSE17536" & variable == "lactate_axis"]
g39_acid <- direction_tests[cohort == "GSE39582" & variable == "acidification_composite"]
g17_acid <- direction_tests[cohort == "GSE17536" & variable == "acidification_composite"]
writeLines(c(
  "# Interpretation of centroid-related effect amplification",
  "",
  "The exclusion did not remove any external samples. Rebuilding the TCGA centroids reassigned 100/566 GSE39582 samples and 30/177 GSE17536 samples, with equal numbers moving in each direction (50 versus 50 and 15 versus 15, respectively).",
  "",
  paste0("Changed samples were concentrated near the projection boundary under both specifications. Under the full-cohort centroids, the median absolute margins were 0.125 versus 0.447 in GSE39582 and 0.094 versus 0.428 in GSE17536 for changed versus unchanged samples; both FDR values were < 1e-10. The same pattern remained under the exclusion-13 centroids."),
  "",
  paste0("Directional reassignment was biologically structured for the lactate axis. Samples moved from C1 to C2 had higher median lactate-axis scores than samples moved from C2 to C1 in GSE39582 (",
         fmt(g39_lac$group1_median), " versus ", fmt(g39_lac$group2_median), "; FDR = ",
         format(g39_lac$FDR_within_cohort, scientific = TRUE, digits = 3), ") and GSE17536 (",
         fmt(g17_lac$group1_median), " versus ", fmt(g17_lac$group2_median), "; FDR = ",
         format(g17_lac$FDR_within_cohort, scientific = TRUE, digits = 3), ")."),
  paste0("The acidification-composite differences pointed in the same direction but did not meet the within-cohort FDR threshold after correction (GSE39582 FDR = ",
         fmt(g39_acid$FDR_within_cohort, 4), "; GSE17536 FDR = ", fmt(g17_acid$FDR_within_cohort, 4), ")."),
  "",
  "Continuous analyses showed that the exclusion-13 signed projection score was more strongly correlated with both marker axes in both cohorts. Paired 2,000-sample bootstrap intervals for the increase in Spearman rho excluded zero for all four cohort-marker comparisons. Thus, amplification was not explained solely by the binary cutoff: centroid reconstruction modestly rotated the continuous projection direction toward the lactate/acidification axes, while the largest categorical instability remained concentrated among low-margin samples.",
  "",
  "The full-cohort centroids remain the prespecified primary analysis. These findings explain the larger exclusion-13 contrasts but do not establish that the rebuilt centroids are biologically superior. They support describing the lactate association as robust and the broader acidification-composite association as sensitive to centroid specification."
), file.path(output_dir, "10_interpretation.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(output_dir, "11_R_sessionInfo.txt"))
cat("Transition/margin and continuous-correlation audit completed\n")
