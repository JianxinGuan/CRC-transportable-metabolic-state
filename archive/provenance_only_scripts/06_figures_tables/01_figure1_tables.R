options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

out <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_results63")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

inp <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_results31_32/inputs")
audit <- fread(file.path(inp, "meta_context_inclusion_audit.tsv"))
intersection <- fread(file.path(inp, "intersection_summary.tsv"))
tcga_summary <- fread(file.path(inp, "tcga_expression_filter_summary.tsv"))
final_summary <- fread(file.path(inp, "final_gene_set_summary.tsv"))
meta_robust <- fread(file.path(inp, "meta_set_robustness_summary.tsv"))
crispr_cond <- fread(file.path(inp, "crispr_condition_summary.tsv"))
map_summary <- fread(file.path(inp, "ensembl_mapping_summary.tsv"))
c174 <- fread(file.path(inp, "dual_evidence_174.tsv"))
g68 <- fread(file.path(inp, "frozen_68.tsv"))

# Formal funnel: TCGA and external CRC filters are parallel branches from the 174 candidates.
funnel <- data.table(
  step_id = c("F1", "F2", "F3a", "F3b", "F4"),
  stage = c(
    "Strict human acidosis meta-analysis",
    "One-to-one human-CRISPR intersection",
    "TCGA expression detectability/variance filter",
    "External CRC cross-platform availability filter",
    "Intersection of TCGA-pass and external-pass genes"
  ),
  input_n = c(NA_integer_, 1500L, 174L, 174L, NA_integer_),
  excluded_n = c(NA_integer_, 1326L, 4L, 104L, NA_integer_),
  retained_n = c(1500L, 174L, 170L, 70L, 68L),
  rule = c(
    "k>=4; random-effects BH FDR<0.05; direction consistency>=0.75; leave-one-dataset-out sign stable",
    "Ensembl release 113 one-to-one ortholog; CRISPR dependency supported in >=2 stress conditions",
    "TPM>1 in >=80% of 647 primary tumors and non-zero variance",
    "Available after the frozen external-CRC cross-platform filtering chain",
    "Set intersection; no survival, CMS, drug-response or projection outcome used for selection"
  ),
  branch = c("discovery", "discovery", "parallel_TCGA", "parallel_external", "final_intersection"),
  outcome_used = "No"
)
fwrite(funnel, file.path(out, "Table_174_to_68_formal_funnel.tsv"), sep = "\t")

# Dataset/context table used by Figure 1 and Table 1 drafting.
context <- audit[, .(
  dataset, contrast, contrast_family, cell_line, analysis_method, platform,
  role = meta_role
)]
context[, formal_primary := role == "primary"]
fwrite(context, file.path(out, "Figure1_GEO_contexts_and_methods.tsv"), sep = "\t")

# Frozen numerical audit.
candidate_stats <- data.table(
  metric = c(
    "primary_GEO_datasets", "primary_experimental_contexts", "RNAseq_contexts", "microarray_contexts",
    "robust_meta_genes", "strict_meta_genes", "strict_one_to_one_dual_evidence",
    "strict_all_mapped_sensitivity", "dual_evidence_up_in_acidosis", "dual_evidence_down_in_acidosis",
    "TCGA_samples", "TCGA_pass_genes", "external_CRC_pass_genes", "frozen_final_genes",
    "final68_subset_of174"
  ),
  value = c(
    uniqueN(audit[meta_role == "primary", dataset]), nrow(audit[meta_role == "primary"]),
    sum(audit[meta_role == "primary", grepl("RNA-seq", platform)]),
    sum(audit[meta_role == "primary", !grepl("RNA-seq", platform)]),
    intersection[set == "robust_meta", n], intersection[set == "strict_meta", n],
    intersection[set == "strict_x_one_to_one", n], intersection[set == "strict_x_all_mapped", n],
    sum(c174$random_effect > 0), sum(c174$random_effect < 0),
    as.numeric(tcga_summary[metric == "tcga_samples", value]),
    as.numeric(tcga_summary[metric == "tcga_pass", value]),
    as.numeric(final_summary$stage17_external), as.numeric(final_summary$final_intersection),
    all(toupper(g68[[1]]) %chin% toupper(c174$gene_symbol))
  ),
  source = c(
    rep("stage16 meta_context_inclusion_audit.tsv", 4),
    rep("stage58 intersection_summary.tsv", 4),
    rep("stage58 dual_evidence_174_candidates_full.tsv", 2),
    rep("stage21 TCGA_expression_filter_summary.tsv", 2),
    rep("stage22 CRC_final_gene_set_summary.tsv", 2),
    "stage58 frozen gene lists"
  )
)
fwrite(candidate_stats, file.path(out, "Results_3.1_3.2_frozen_number_audit.tsv"), sep = "\t")

# Figure panel data.
meta_counts <- data.table(
  set = factor(c("Robust meta", "Strict meta", "Strict x 1:1 ortholog", "Frozen CRC set"),
               levels = rev(c("Robust meta", "Strict meta", "Strict x 1:1 ortholog", "Frozen CRC set"))),
  n = c(1647, 1500, 174, 68),
  class = c("Human meta", "Human meta", "Cross-species", "CRC filtered")
)

cond_labels <- c(
  low_glucose_aa = "Low glucose +\nlow NEAA",
  low_glucose = "Low glucose",
  hypoxia = "Hypoxia",
  lactic_acid = "Lactic acid",
  low_NEAA = "Low NEAA"
)
crispr_plot <- copy(crispr_cond)
crispr_plot[, label := factor(cond_labels[condition], levels = rev(unname(cond_labels)))]

theme_fig <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5), axis.title = element_text(size = 9),
        legend.position = "bottom", plot.margin = margin(7, 7, 7, 7))

pA <- ggplot() +
  annotate("rect", xmin = 0.07, xmax = 0.93,
           ymin = c(0.665, 0.365, 0.065), ymax = c(0.895, 0.595, 0.295),
           fill = c("#D8E8F2", "#E6DFF0", "#DCECDD"),
           color = "#3F4A52", linewidth = 0.75) +
  annotate("text", x = 0.50, y = c(0.780, 0.480, 0.180),
           label = c("Human acidosis transcriptomes\n7 datasets / 11 contrasts",
                     "Metabolic-stress CRISPR\n5 stress conditions",
                     "CRC-specific expression filtering\nTCGA / external cohorts"),
           color = "#20262B", size = 2.85, lineheight = 1.10) +
  annotate("segment", x = 0.50, xend = 0.50,
           y = c(0.645, 0.345), yend = c(0.615, 0.315),
           arrow = arrow(length = unit(0.11, "in"), type = "closed"),
           linewidth = 0.72, color = "#3F4A52") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.95), clip = "off") +
  labs(title = "A  Cross-species discovery framework") + theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0), plot.margin = margin(7, 7, 7, 7))

ctx <- audit[meta_role == "primary", .N, by = .(dataset, analysis_method)]
pB <- ggplot(ctx, aes(x = reorder(dataset, N), y = N, fill = analysis_method)) +
  geom_col(width = 0.72) + coord_flip() +
  geom_text(aes(label = N), hjust = -0.2, size = 3) +
  scale_fill_manual(values = c(DESeq2 = "#4477AA", limma = "#CC6677")) +
  scale_y_continuous(expand = expansion(mult = c(0, .15)), breaks = 0:4) +
  labs(title = "B  Human discovery contexts", subtitle = "Independent contrasts; matrices were not pooled", x = NULL, y = "Contrasts", fill = "Method") + theme_fig

pC <- ggplot(meta_counts, aes(x = set, y = n, fill = class)) +
  geom_col(width = 0.7) + coord_flip() + geom_text(aes(label = n), hjust = -0.12, size = 3.2) +
  scale_fill_manual(values = c("Human meta" = "#4477AA", "Cross-species" = "#AA4499", "CRC filtered" = "#228833")) +
  scale_y_continuous(expand = expansion(mult = c(0, .16))) +
  labs(title = "C  Candidate reduction", subtitle = "Random-effects meta-analysis and functional constraint", x = NULL, y = "Genes", fill = NULL) + theme_fig

pD <- ggplot(crispr_plot, aes(x = label, y = FDR05_negative)) +
  geom_col(fill = "#AA4499", width = 0.7) + coord_flip() +
  geom_text(aes(label = FDR05_negative), hjust = -0.12, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, .16))) +
  labs(title = "D  CRISPR dependency hits", subtitle = "Negative LFC; FDR < 0.05; 18,202 tested", x = NULL, y = "Genes") + theme_fig

fdat <- data.table(
  x = c(1, 2, 3, 3, 4), y = c(3, 3, 4, 2, 3),
  label = c("Strict human Meta\n1,500", "1:1 CRISPR ortholog\n174", "TCGA pass\n170", "External CRC pass\n70", "Frozen intersection\n68"),
  fill = c("#DDEAF3", "#E7E1F2", "#E2EFE3", "#F6E4D7", "#CFE8D5")
)
pE <- ggplot(fdat, aes(x, y)) +
  geom_segment(data = data.table(x = c(1.25, 2.25, 2.25, 3.25, 3.25), y = c(3, 3, 3, 4, 2),
                                 xend = c(1.75, 2.75, 2.75, 3.75, 3.75), yend = c(3, 4, 2, 3, 3)),
               aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE,
               arrow = arrow(length = unit(.10, "in")), linewidth = .45, color = "#56616A") +
  geom_label(aes(label = label, fill = fill), size = 3.05, label.size = .25, label.r = unit(.08, "lines"), lineheight = 1.0) +
  scale_fill_identity() + coord_cartesian(xlim = c(.65, 4.35), ylim = c(1.3, 4.7), clip = "off") +
  labs(title = "E  Formal 174-to-68 funnel", subtitle = "TCGA and external CRC filters were parallel; final set is their intersection") +
  theme_void(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 8.5), plot.margin = margin(7,7,7,7))

fig <- (pA | pB) / (pC | pD) / pE + plot_layout(heights = c(.8, 1, .9))
ggsave(file.path(out, "Figure1_analysis_version.pdf"), fig, width = 180, height = 230, units = "mm", device = cairo_pdf)
ggsave(file.path(out, "Figure1_analysis_version_300dpi.png"), fig, width = 180, height = 230, units = "mm", dpi = 300, bg = "white")

writeLines(c(
  "Figure 1. Cross-species discovery and CRC-specific freezing of the 68-gene signature.",
  "(A) Study design integrating independently analyzed human acidosis transcriptomic contexts, five metabolic-stress CRISPR conditions, and CRC expression filters.",
  "(B) Seven primary GEO datasets contributed 11 experimental contrasts; RNA-seq and microarray contrasts were analyzed with DESeq2 and limma, respectively, before meta-analysis.",
  "(C) Random-effects integration yielded 1,647 robust and 1,500 strict human candidates. Intersection with Ensembl release 113 one-to-one CRISPR orthologs supported in at least two stress conditions produced 174 dual-evidence candidates, of which 68 entered the frozen CRC signature.",
  "(D) Numbers of negative-dependency hits (negative MAGeCK LFC, FDR<0.05) in each CRISPR stress condition. Bars describe the original condition-level screen and are not additive because genes can be supported by multiple conditions.",
  "(E) Formal filtering logic. TCGA expression filtering retained 170 of 174 genes, whereas the independently defined external-CRC availability filter retained 70; their intersection contained 68 genes. Survival, CMS, drug-response and external-projection outcomes were not used for gene selection.",
  "Abbreviations: CRC, colorectal cancer; FDR, false discovery rate; LFC, log fold change; NEAA, non-essential amino acids; TCGA, The Cancer Genome Atlas."
), file.path(out, "Figure1_legend_draft.md"))

capture.output(sessionInfo(), file = file.path(out, "R_sessionInfo.txt"))
cat("Figure 1 and tables completed\n")
