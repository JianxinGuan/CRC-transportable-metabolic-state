options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(20260830)
inp <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_results33_34/inputs")
out <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_results64")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

stab <- fread(file.path(inp, "stability.tsv"))
lab <- fread(file.path(inp, "k2_labels.tsv"))
score <- fread(file.path(inp, "state_scores.tsv"))
cm0 <- fread(file.path(inp, "k2_consensus.tsv.gz"))
proj_s <- fread(file.path(inp, "projection_samples.tsv.gz"))
proj <- fread(file.path(inp, "projection_summary.tsv"))
sil <- fread(file.path(inp, "silhouette_comparison.tsv"))
sil_ci <- fread(file.path(inp, "silhouette_CI.tsv"))
new_s <- fread(file.path(inp, "new_projection_samples.tsv.gz"))
new_proj <- fread(file.path(inp, "new_projection_summary.tsv"))
margin <- fread(file.path(inp, "margin_comparison.tsv"))

# Freeze direction and sample-level table.
tcga <- merge(lab, score[, .(sample_id, acidosis_signed_score)], by = "sample_id")
tcga[, state_end := factor(paste0("C", cluster), levels = c("C1", "C2"))]
setorder(tcga, acidosis_signed_score)
tcga[, score_rank := seq_len(.N)]
fwrite(tcga, file.path(out, "Figure2_TCGA_continuous_state_scores.tsv"), sep = "\t")

# Consensus matrix ordered by K=2 labels and within-cluster hierarchical order.
ids <- cm0[[1]]
cm <- as.matrix(cm0[, -1])
rownames(cm) <- ids
storage.mode(cm) <- "double"
common <- intersect(lab$sample_id, rownames(cm))
ord <- unlist(lapply(1:2, function(cl) {
  z <- lab[cluster == cl & sample_id %chin% common, sample_id]
  if (length(z) > 2) z[hclust(as.dist(1 - cm[z, z, drop = FALSE]), method = "average")$order] else z
}))
cm <- cm[ord, ord, drop = FALSE]
cm_long <- data.table(
  row = rep(seq_len(nrow(cm)), times = ncol(cm)),
  col = rep(seq_len(ncol(cm)), each = nrow(cm)),
  consensus = as.vector(cm)
)
boundary <- sum(lab$cluster[match(ord, lab$sample_id)] == 1)

# Stability summary data.
stab[, smallest_cluster := vapply(strsplit(cluster_sizes, ";", fixed = TRUE), function(x) min(as.integer(x)), integer(1))]
stab[, selection_score := mean_consensus - sd_consensus]
stab[, selected := k == 2]
fwrite(stab, file.path(out, "Figure2_K2_K6_resampling_summary.tsv"), sep = "\t")

# Frozen silhouette table with bootstrap confidence intervals.
sil_main <- sil[analysis %chin% c("Frozen68_K2", "Deletion47_K2", "E2F_G2M_full_K2", "E2F_G2M_exclusive_K2")]
sil_main <- merge(sil_main, sil_ci, by = "analysis", all.x = TRUE)
sil_main[, label := factor(analysis,
  levels = c("Frozen68_K2", "Deletion47_K2", "E2F_G2M_full_K2", "E2F_G2M_exclusive_K2"),
  labels = c("Frozen 68", "Deletion 47", "E2F/G2M full", "E2F/G2M exclusive"))]
fwrite(sil_main, file.path(out, "Figure2_TCGA_silhouette_comparison.tsv"), sep = "\t")

# External projection summaries.
ext <- proj[cohort != "TCGA_COADREAD"]
ext_long <- melt(ext[, .(cohort, cluster1_n, cluster2_n)], id.vars = "cohort",
                 variable.name = "state_end", value.name = "n")
ext_long[, state_end := factor(state_end, levels = c("cluster1_n", "cluster2_n"), labels = c("C1", "C2"))]
ext_long[, percent := 100 * n / sum(n), by = cohort]
fwrite(ext_long, file.path(out, "Figure2_external_projection_subtype_sizes.tsv"), sep = "\t")

sil_ext <- rbind(
  ext[, .(cohort, centroid = "Original TCGA centroids", silhouette_mean)],
  new_proj[, .(cohort, centroid = "Exclusion-13 centroids", silhouette_mean = new_silhouette_mean)]
)
sil_ext[, centroid := factor(centroid, levels = c("Original TCGA centroids", "Exclusion-13 centroids"))]
fwrite(sil_ext, file.path(out, "Figure2_external_silhouette_sensitivity.tsv"), sep = "\t")

margin_samples <- new_s[, .(cohort, changed = factor(ifelse(changed, "Changed", "Unchanged"),
                                                      levels = c("Changed", "Unchanged")), old_margin)]
fwrite(margin_samples, file.path(out, "Figure2_projection_margin_sample_level.tsv.gz"), sep = "\t")

# Complete frozen-number audit.
audit <- data.table(
  metric = c("TCGA_n", "TCGA_C1_n", "TCGA_C2_n", "TCGA_frozen68_silhouette",
             "GSE39582_n", "GSE39582_C1_n", "GSE39582_C2_n", "GSE39582_original_silhouette",
             "GSE17536_n", "GSE17536_C1_n", "GSE17536_C2_n", "GSE17536_original_silhouette",
             "GSE39582_new_silhouette", "GSE17536_new_silhouette",
             "GSE39582_changed_n", "GSE17536_changed_n", "GSE39582_margin_FDR", "GSE17536_margin_FDR"),
  value = c(proj[cohort == "TCGA_COADREAD", n], proj[cohort == "TCGA_COADREAD", cluster1_n],
            proj[cohort == "TCGA_COADREAD", cluster2_n], proj[cohort == "TCGA_COADREAD", silhouette_mean],
            ext[cohort == "GSE39582", n], ext[cohort == "GSE39582", cluster1_n], ext[cohort == "GSE39582", cluster2_n], ext[cohort == "GSE39582", silhouette_mean],
            ext[cohort == "GSE17536", n], ext[cohort == "GSE17536", cluster1_n], ext[cohort == "GSE17536", cluster2_n], ext[cohort == "GSE17536", silhouette_mean],
            new_proj[cohort == "GSE39582", new_silhouette_mean], new_proj[cohort == "GSE17536", new_silhouette_mean],
            new_proj[cohort == "GSE39582", changed_n], new_proj[cohort == "GSE17536", changed_n],
            margin[cohort == "GSE39582", FDR], margin[cohort == "GSE17536", FDR]),
  interpretation = c("Discovery expression cohort", "Original K=2", "Original K=2", "Mean silhouette",
                     "Expression/projection cohort", "Original centroid projection", "Original centroid projection", "Mean silhouette",
                     "Expression/projection cohort", "Original centroid projection", "Original centroid projection", "Mean silhouette",
                     "Mean silhouette after rebuilding centroids", "Mean silhouette after rebuilding centroids",
                     "Original vs exclusion-13 centroid label changes", "Original vs exclusion-13 centroid label changes",
                     "BH-adjusted Wilcoxon test", "BH-adjusted Wilcoxon test")
)
fwrite(audit, file.path(out, "Results_3.3_3.4_frozen_number_audit.tsv"), sep = "\t")

# Plot styling.
cols <- c(C1 = "#2678A5", C2 = "#D65F5F")
theme_fig <- theme_minimal(base_size = 9.5) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8), axis.title = element_text(size = 8.5),
        legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(6, 6, 6, 6))

pA <- ggplot(stab, aes(k, selection_score)) +
  geom_line(color = "#69747C", linewidth = .65) +
  geom_point(aes(fill = selected), shape = 21, size = 3, color = "#35424B") +
  geom_text(aes(label = paste0("min n=", smallest_cluster)), vjust = -1.0, size = 2.5) +
  scale_fill_manual(values = c(`TRUE` = "#2A9D8F", `FALSE` = "white"), guide = "none") +
  scale_x_continuous(breaks = 2:6) +
  scale_y_continuous(expand = expansion(mult = c(.08, .20))) +
  labs(title = "A  K=2-6 resampling summary", subtitle = "Selection score: mean consensus - SD", x = "K", y = "Selection score") + theme_fig

pB <- ggplot(cm_long, aes(col, row, fill = consensus)) +
  geom_raster() +
  geom_vline(xintercept = boundary + .5, color = "white", linewidth = .45) +
  geom_hline(yintercept = boundary + .5, color = "white", linewidth = .45) +
  scale_fill_gradientn(colors = c("#F6F7F8", "#B9D5E5", "#2678A5"), limits = c(0, 1)) +
  scale_y_reverse() +
  coord_fixed(expand = FALSE) +
  labs(title = "B  K=2 consensus", subtitle = "TCGA: C1 n=214; C2 n=433", x = NULL, y = NULL) +
  theme_void(base_size = 9.5) + theme(plot.title = element_text(face = "bold", size = 10.5),
                                     plot.subtitle = element_text(size = 8), legend.position = "none",
                                     plot.margin = margin(6, 6, 6, 6))

pC <- ggplot(tcga, aes(score_rank, acidosis_signed_score, color = state_end)) +
  geom_hline(yintercept = 0, linetype = 2, color = "#89939A", linewidth = .35) +
  geom_point(size = 1.05, alpha = .8) +
  scale_color_manual(values = cols) +
  labs(title = "C  Continuous 68-gene state score", subtitle = "K=2 labels mark the two ends of a continuous ranking",
       x = "Samples ordered by score", y = "Signed state score", color = NULL) + theme_fig

pD <- ggplot(sil_main, aes(label, silhouette_mean, fill = label)) +
  geom_hline(yintercept = .1, linetype = 2, color = "#666666", linewidth = .4) +
  geom_col(width = .68) +
  geom_errorbar(aes(ymin = CI95_low, ymax = CI95_high), width = .16, linewidth = .45) +
  geom_text(aes(label = sprintf("%.3f", silhouette_mean)), vjust = -1.1, size = 2.6) +
  scale_fill_manual(values = c("#2678A5", "#6AA6C8", "#9D79BC", "#76518F"), guide = "none") +
  scale_y_continuous(limits = c(0, .115), expand = expansion(mult = c(0, .03))) +
  labs(title = "D  TCGA geometric separation", subtitle = "Mean silhouette; 2,000-sample bootstrap 95% CI", x = NULL, y = "Mean silhouette") +
  theme_fig + theme(axis.text.x = element_text(angle = 25, hjust = 1))

pE <- ggplot(ext_long, aes(cohort, percent, fill = state_end)) +
  geom_col(width = .65) +
  geom_text(aes(label = paste0(n, " (", sprintf("%.1f", percent), "%)")),
            position = position_stack(vjust = .5), color = "white", size = 2.7) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  labs(title = "E  External projected state sizes", subtitle = "Frozen original TCGA centroids", x = NULL, y = "Samples", fill = NULL) + theme_fig

pF <- ggplot(sil_ext, aes(cohort, silhouette_mean, fill = centroid)) +
  geom_col(position = position_dodge(.72), width = .65) +
  geom_text(aes(label = sprintf("%.3f", silhouette_mean)), position = position_dodge(.72), vjust = -0.7, size = 2.6) +
  scale_fill_manual(values = c("#6E8FA8", "#2A9D8F")) +
  scale_y_continuous(limits = c(0, .15), expand = expansion(mult = c(0, .03))) +
  labs(title = "F  External silhouette sensitivity", subtitle = "Original vs exclusion-13-derived centroids", x = NULL, y = "Mean silhouette") + theme_fig

pG <- ggplot(margin_samples, aes(changed, old_margin, fill = changed)) +
  geom_boxplot(width = .62, outlier.shape = NA, linewidth = .4) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_fill_manual(values = c(Changed = "#D65F5F", Unchanged = "#7A9E7E"), guide = "none") +
  coord_cartesian(ylim = c(0, quantile(margin_samples$old_margin, .97))) +
  labs(title = "G  Label changes concentrate at low margins", subtitle = "Old margin = |correlation to C1 - correlation to C2|",
       x = NULL, y = "Original projection margin") + theme_fig + theme(axis.text.x = element_text(angle = 20, hjust = 1))

fig <- (pA | pB) / (pC | pD) / (pE | pF) / pG + plot_layout(heights = c(1, 1, .95, 1.05))
ggsave(file.path(out, "Figure2_analysis_version.pdf"), fig, width = 180, height = 285, units = "mm", device = cairo_pdf)
ggsave(file.path(out, "Figure2_analysis_version_300dpi.png"), fig, width = 180, height = 285, units = "mm", dpi = 300, bg = "white")

writeLines(c(
  "Figure 2. The frozen 68-gene signature delineates a continuous TCGA expression state that is transferable to independent CRC cohorts.",
  "(A) Summary of K=2-6 resampling results. K=2 maximized the predefined mean-consensus-minus-standard-deviation score and avoided the very small clusters observed for K>=3. The metric is a resampling summary rather than evidence of discrete biological classes.",
  "(B) K=2 consensus matrix for 647 TCGA-COAD/READ primary tumors ordered by final labels (C1, n=214; C2, n=433).",
  "(C) Signed 68-gene state scores across TCGA samples ordered continuously. K=2 labels are shown as operational labels for the two ends of the ranking.",
  "(D) Mean TCGA silhouette widths for the frozen 68-gene set, the 47-gene deletion sensitivity set and two pure-proliferation controls. Error bars are percentile 95% confidence intervals from 2,000 bootstrap replicates; the dashed line marks 0.1. Silhouette describes geometric separation, not classification accuracy.",
  "(E) C1/C2 sizes after projecting the frozen original TCGA centroids into GSE39582 and GSE17536.",
  "(F) Mean external silhouette widths obtained with the original TCGA centroids and with centroids rebuilt after excluding 13 source/vial-associated TCGA samples.",
  "(G) Original projection margins among samples whose labels changed or remained unchanged after centroid rebuilding. Label changes were concentrated among low-margin samples (BH FDR=2.21e-24 for GSE39582 and 3.75e-11 for GSE17536; two-sided Wilcoxon rank-sum tests).",
  "Abbreviations: BH, Benjamini-Hochberg; CRC, colorectal cancer; TCGA, The Cancer Genome Atlas."
), file.path(out, "Figure2_legend_draft.md"))

capture.output(sessionInfo(), file = file.path(out, "R_sessionInfo.txt"))
cat("Figure 2 completed\n")
