options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")) else getwd()
target <- if (length(args) >= 1) args[1] else script_dir

overlap_plot <- fread(file = file.path(target, "Figure4_marker_overlap_audit.tsv"))
k2 <- fread(file = file.path(target, "Figure4_TCGA_K2_marker_differences.tsv"))
cors <- fread(file = file.path(target, "Figure4_TCGA_unadjusted_correlations.tsv"))
forest <- fread(file = file.path(target, "Figure4_adjusted_models.tsv"))
extplot <- fread(file = file.path(target, "Figure4_external_marker_effects_two_centroids_source.tsv"))

marker_order_7 <- rev(c("Carbonic anhydrases", "Acid/proton transport", "Lactate axis",
                        "CA9/CA12", "Acidification composite", "CA9", "CA12"))
marker_order_5 <- rev(c("Carbonic anhydrases", "Acid/proton transport", "Lactate axis",
                        "CA9/CA12", "Acidification composite"))
overlap_order <- rev(c("Carbonic anhydrases", "Acid/proton transport", "Lactate axis", "CA9/CA12"))

overlap_plot[, label := factor(as.character(label), levels = overlap_order)]
k2[, label := factor(as.character(label), levels = marker_order_7)]
cors[, label := factor(as.character(label), levels = marker_order_7)]
forest[, label := factor(as.character(label), levels = marker_order_5)]
forest[, analysis := factor(as.character(analysis),
                            levels = c("All evaluable", "Exclude 13"),
                            labels = c("Full TCGA cohort", "Exclusion-13 sensitivity set"))]

extplot <- extplot[outcome %chin% c("acidification_composite", "lactate_axis")]
extplot[, `:=`(
  label = factor(ifelse(outcome == "lactate_axis", "Lactate axis", "Acidification composite"),
                 levels = c("Acidification composite", "Lactate axis")),
  cohort = factor(cohort, levels = c("GSE39582", "GSE17536")),
  centroid = factor(centroid_specification,
                    levels = c("full_cohort_fixed_centroids", "exclusion13_rebuilt_centroids"),
                    labels = c("Full-cohort centroids", "Exclusion-13 centroids")),
  effect = delta_C2_minus_C1,
  FDR = FDR_within_cohort_and_specification,
  significant = FDR_within_cohort_and_specification < 0.05,
  point_fill = ifelse(centroid_specification == "exclusion13_rebuilt_centroids", "white",
                      ifelse(cohort == "GSE39582", "#2A9D8F", "#E07A5F"))
)]

# This is the formal source table for Figure 4E. It intentionally retains both
# the prespecified full-cohort analysis and the exclusion-13 sensitivity analysis.
fwrite(extplot[, .(cohort, centroid_specification, centroid, outcome, label, n,
                   cluster1_n, cluster2_n, genes_detected, C1_mean, C2_mean,
                   effect_C2_minus_C1 = effect, p_value, FDR, significant)],
       file.path(target, "Figure4_external_marker_effects.tsv"), sep = "\t")

theme_fig <- theme_minimal(base_size = 9.5) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8),
        axis.title = element_text(size = 8.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.margin = margin(6, 6, 6, 6))

p4A <- ggplot(overlap_plot, aes(n_detected, label)) +
  geom_col(fill = "#5B8E7D", width = .65) +
  geom_text(aes(label = paste0(n_detected, " detected\n", n_overlap68, " overlap")),
            hjust = -.08, size = 2.8) +
  scale_x_continuous(limits = c(0, 18), expand = expansion(mult = c(0, .02))) +
  labs(title = "A  Independent marker-set audit",
       subtitle = "Four constituent sets detected; zero overlap with the frozen 68 genes",
       x = "Genes", y = NULL) + theme_fig

p4B <- ggplot(k2, aes(delta_C2_minus_C1, label, fill = FDR < .05)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#777777") +
  geom_col(width = .65) +
  scale_fill_manual(values = c(`TRUE` = "#C95D63", `FALSE` = "#B9B9B9"), guide = "none") +
  labs(title = "B  TCGA categorical differences", subtitle = "Unadjusted effect = C2 - C1",
       x = "Marker-score difference", y = NULL) + theme_fig

p4C <- ggplot(cors, aes(rho, label, fill = FDR < .05)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#777777") +
  geom_col(width = .65) +
  scale_fill_manual(values = c(`TRUE` = "#557A95", `FALSE` = "#B9B9B9"), guide = "none") +
  labs(title = "C  TCGA unadjusted correlations",
       subtitle = "Spearman correlation with continuous state score",
       x = "Spearman rho", y = NULL) + theme_fig

p4D <- ggplot(forest, aes(beta, label, color = analysis)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#777777") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = .18,
                 position = position_dodge(.45)) +
  geom_point(size = 2.3, position = position_dodge(.45)) +
  scale_color_manual(values = c("#2F6690", "#D77A61")) +
  labs(title = "D  Fully adjusted TCGA models",
       subtitle = "Beta per 1-unit state-score increase",
       x = "Adjusted beta (95% CI)", y = NULL) + theme_fig +
  theme(legend.text = element_text(size = 7.4))

# Paired point-range display makes the prespecified and sensitivity estimates
# visible without implying that only the significant sensitivity result is primary.
p4E <- ggplot(extplot, aes(effect, cohort, group = cohort, color = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#777777") +
  geom_line(aes(group = cohort), color = "#A6A6A6", linewidth = .55) +
  geom_point(aes(shape = centroid, fill = point_fill), size = 3.0, stroke = .8,
             show.legend = TRUE) +
  geom_text(aes(x = effect + .012, label = ifelse(significant, "*", "ns")), size = 2.7,
            hjust = 0, show.legend = FALSE) +
  facet_wrap(~label, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c(GSE39582 = "#2A9D8F", GSE17536 = "#E07A5F"), guide = "none") +
  scale_fill_identity(guide = "none") +
  scale_shape_manual(values = c("Full-cohort centroids" = 16, "Exclusion-13 centroids" = 21)) +
  scale_x_continuous(limits = c(0, .49), expand = expansion(mult = c(.01, .03))) +
  labs(title = "E  External centroid sensitivity",
       subtitle = "C2 - C1; filled = full-cohort; open = exclusion-13; * FDR < 0.05",
       x = "Marker-score difference", y = NULL, shape = NULL) + theme_fig +
  theme(strip.text = element_text(face = "bold", size = 7.8),
        legend.text = element_text(size = 7.2), legend.key.width = unit(3.2, "mm"))

p4F <- ggplot() +
  annotate("rect", xmin = .04, xmax = .96, ymin = .08, ymax = .92,
           fill = "#F3F5F6", color = "#64717A", linewidth = .6) +
  annotate("text", x = .5, y = .73,
           label = "Patient-side transcriptomic\nassociation",
           fontface = "bold", size = 3.4, lineheight = .95) +
  annotate("text", x = .5, y = .48,
           label = "Centroid-robust lactate signal;\nbroader acidification association\nis specification-sensitive",
           size = 2.75, lineheight = .95) +
  annotate("text", x = .5, y = .20,
           label = "Does not directly measure tumor pH\nor establish causal adaptation",
           size = 2.7, lineheight = .95, color = "#A23B3B") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "F  Interpretation boundary") + theme_void() +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.margin = margin(6, 6, 6, 6))

fig4 <- (p4A | p4B) / (p4C | p4D) / (p4E | p4F) +
  plot_layout(heights = c(1, 1.05, 1.08))

ggsave(file.path(target, "Figure4_analysis_version.pdf"), fig4,
       width = 200, height = 235, units = "mm", device = cairo_pdf)
ggsave(file.path(target, "Figure4_analysis_version_300dpi.png"), fig4,
       width = 200, height = 235, units = "mm", dpi = 300, bg = "white")

legend_text <- c(
  "Figure 4. Independent acid-handling and lactate-related transcriptomic markers provide qualified patient-side support for an acidosis-linked interpretation.",
  "(A) Four constituent marker collections were detected and had zero overlap with the frozen 68-gene feature. Together with the derived acidification composite, these constituted five prespecified outcomes.",
  "(B) Unadjusted TCGA C2-minus-C1 marker-score differences.",
  "(C) Unadjusted Spearman correlations between marker scores and the continuous state score.",
  "(D) Fully adjusted TCGA linear models in the full evaluable cohort and after excluding 13 source/vial-associated samples. Models adjusted for proliferation, Hallmark hypoxia, an xCell-derived microenvironment proxy, CMS, COAD/READ, and stage.",
  "(E) External C2-minus-C1 differences under the prespecified full-cohort TCGA centroids and exclusion-13 rebuilt centroids. Lactate-axis differences remained significant in both cohorts under both specifications. The acidification-composite direction was consistent, but significance was confined to the exclusion-13 specification. Benjamini-Hochberg correction was applied across five prespecified marker outcomes within each cohort and centroid specification.",
  "(F) Interpretation boundary. The results support a centroid-robust lactate-related transcriptional association and a specification-sensitive broader acidification association; they do not directly measure intratumoral pH or establish causal adaptation to acidosis."
)
writeLines(legend_text, file.path(target, "Figure4_legend_draft.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(target, "R_sessionInfo_Figure4_rebuild.txt"))
cat("Updated Figure 4 completed\n")
