suppressPackageStartupMessages({library(data.table); library(cluster); library(survival); library(GSVA)})
options(stringsAsFactors = FALSE)
out <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage61/10_external_projection_exclusion13_centroids")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(p) fread(p, sep = "\t", data.table = TRUE)
read_gz <- function(p) as.data.table(read.delim(gzfile(p), check.names = FALSE))
z_rows <- function(m) { z <- t(scale(t(m))); z[!is.finite(z)] <- 0; z }
ari <- function(x, y) {
  tab <- table(x, y); n <- sum(tab); c2 <- function(z) z * (z - 1) / 2
  a <- sum(c2(tab)); b <- sum(c2(rowSums(tab))); c <- sum(c2(colSums(tab)))
  ex <- b * c / c2(n); den <- (b + c) / 2 - ex
  if (den == 0) NA_real_ else (a - ex) / den
}
nmi <- function(x, y) {
  tab <- table(x, y); p <- tab / sum(tab); px <- rowSums(p); py <- colSums(p); e <- outer(px, py); z <- p > 0
  mi <- sum(p[z] * log(p[z] / e[z])); mi / sqrt((-sum(px[px > 0] * log(px[px > 0]))) * (-sum(py[py > 0] * log(py[py > 0]))))
}
save_plot <- function(stem, draw, w = 8, h = 6) {
  pdf(paste0(stem, ".pdf"), width = w, height = h, useDingbats = FALSE); draw(); dev.off()
  png(paste0(stem, ".png"), width = w, height = h, units = "in", res = 220); draw(); dev.off()
}

# New TCGA centroids from the 634-sample exclusion rerun.
tc <- read_tsv(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "CRC68_TCGA_gene_zscore.tsv"))
genes <- sub("^\\ufeff", "", tc[[1]])
tm <- as.matrix(tc[, -1]); rownames(tm) <- genes; storage.mode(tm) <- "double"
new_lab <- read_tsv(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage61/07_small13_exclusion_sensitivity/exclusion13_rerun_68gene_K2_labels.tsv"))
ids <- intersect(new_lab$sample_id, colnames(tm)); new_lab <- new_lab[match(ids, sample_id)]; tm <- tm[, ids, drop = FALSE]
cent <- sapply(1:2, function(k) rowMeans(tm[, new_lab$exclusion13_rerun_K2 == k, drop = FALSE]))
colnames(cent) <- c("C1_exclusion13", "C2_exclusion13")
fwrite(data.table(gene = rownames(cent), C1_exclusion13 = cent[, 1], C2_exclusion13 = cent[, 2]),
       file.path(out, "TCGA_exclusion13_centroids.tsv"), sep = "\t")

g395x <- readRDS(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "GSE39582_expr_surv.rds"))
g175x <- readRDS(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "GSE17536_expr_surv.rds"))
old395 <- read_tsv(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "G", "2026_08_25_Projects/4.第一阶段/46_GSE39582_fixed_K2_projection/GSE39582_K2_projected_labels.tsv"))
old175 <- read_tsv(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "G", "2026_08_25_Projects/4.第一阶段/48_GSE17536_final_external_validation/GSE17536_K2_projected_labels.tsv"))
genes395 <- read_tsv(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "G", "2026_08_25_Projects/4.第一阶段/46_GSE39582_fixed_K2_projection/genes_used.tsv"))[[1]]

project_one <- function(x, old, avail, cohort) {
  avail <- intersect(avail, intersect(names(x), rownames(cent)))
  m <- t(as.matrix(x[, avail, drop = FALSE])); sid <- as.character(x$sample); colnames(m) <- sid
  ez <- z_rows(m)
  c1 <- apply(ez, 2, cor, cent[avail, 1]); c2 <- apply(ez, 2, cor, cent[avail, 2])
  pred <- ifelse(c2 > c1, 2L, 1L)
  old_cl <- old$cluster[match(sid, old$sample_id)]
  changed <- pred != old_cl
  dd <- as.dist(1 - cor(ez)); sw <- silhouette(pred, dd)[, "sil_width"]
  ext_cent <- sapply(1:2, function(k) rowMeans(ez[, pred == k, drop = FALSE]))
  centroid_dist <- c(1 - cor(ext_cent[, 1], cent[avail, 1]), 1 - cor(ext_cent[, 2], cent[avail, 2]))
  os_status <- as.numeric(as.character(x$OS.status)); os_time <- as.numeric(x$OS.time)
  ok <- is.finite(os_status) & is.finite(os_time) & os_time > 0
  sv <- survdiff(Surv(os_time[ok], os_status[ok]) ~ pred[ok]); osp <- pchisq(sv$chisq, 1, lower.tail = FALSE)
  sample <- data.table(cohort = cohort, sample_id = sid, old_cluster = old_cl, new_cluster = pred,
                       changed = changed, old_margin = old$margin[match(sid, old$sample_id)],
                       new_cor_C1 = c1, new_cor_C2 = c2, new_margin = abs(c1 - c2), silhouette = sw,
                       OS_status = os_status, OS_time = os_time)
  summary <- data.table(cohort = cohort, n = length(sid), genes = length(avail),
                        old_cluster1_n = sum(old_cl == 1), old_cluster2_n = sum(old_cl == 2),
                        new_cluster1_n = sum(pred == 1), new_cluster2_n = sum(pred == 2),
                        changed_n = sum(changed), agreement = mean(!changed), ARI = ari(pred, old_cl), NMI = nmi(pred, old_cl),
                        new_silhouette_mean = mean(sw), new_silhouette_median = median(sw),
                        new_negative_silhouette_percent = 100 * mean(sw < 0),
                        new_margin_median = median(abs(c1 - c2)),
                        centroid_distance_C1 = centroid_dist[1], centroid_distance_C2 = centroid_dist[2],
                        old_OS_logrank_p = if (cohort == "GSE39582") 0.243480311618268 else 0.748504444798643,
                        new_OS_logrank_p = osp)
  list(sample = sample, summary = summary, z = ez, labels = data.table(sample_id = sid, cluster = pred))
}

r395 <- project_one(g395x, old395, genes395, "GSE39582")
r175 <- project_one(g175x, old175, genes, "GSE17536")
samples <- rbindlist(list(r395$sample, r175$sample))
summary <- rbindlist(list(r395$summary, r175$summary))
fwrite(samples, file.path(out, "external_reprojection_sample_level.tsv.gz"), sep = "\t")
fwrite(summary, file.path(out, "external_reprojection_quality_summary.tsv"), sep = "\t")
fwrite(r395$labels, file.path(out, "GSE39582_exclusion13_centroid_labels.tsv"), sep = "\t")
fwrite(r175$labels, file.path(out, "GSE17536_exclusion13_centroid_labels.tsv"), sep = "\t")

# Changed samples should be concentrated near the old projection boundary.
margin_cmp <- rbindlist(lapply(unique(samples$cohort), function(co) {
  d <- samples[cohort == co]; w <- wilcox.test(d$old_margin[d$changed], d$old_margin[!d$changed], exact = FALSE)
  data.table(cohort = co, changed_n = sum(d$changed), unchanged_n = sum(!d$changed),
             changed_old_margin_median = median(d$old_margin[d$changed]),
             unchanged_old_margin_median = median(d$old_margin[!d$changed]), p_value = w$p.value)
}))
margin_cmp[, FDR := p.adjust(p_value, "BH")]
fwrite(margin_cmp, file.path(out, "changed_sample_margin_comparison.tsv"), sep = "\t")

# Hallmark five-pathway and broader metabolic-direction sensitivity.
g395dt <- read_gz(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "G", "2026_08_25_Projects/4.第一阶段/47_GSE39582_Hallmark_GSVA_validation/GSE39582_Hallmark_GSVA_scores.tsv.gz"))
g395sc <- as.matrix(g395dt[, -1]); rownames(g395sc) <- g395dt[[1]]; storage.mode(g395sc) <- "double"
gmt <- strsplit(readLines(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "hallmark_v2026.1.gmt"), warn = FALSE), "\t")
sets <- setNames(lapply(gmt, function(z) z[-c(1, 2)]), vapply(gmt, `[`, "", 1))
g175genes <- setdiff(names(g175x), c("sample", "OS.status", "OS.time"))
g175m <- t(as.matrix(g175x[, g175genes, drop = FALSE])); colnames(g175m) <- as.character(g175x$sample)
sets175 <- lapply(sets, intersect, rownames(g175m)); sets175 <- sets175[lengths(sets175) >= 10]
g175sc <- gsva(gsvaParam(g175m, sets175, kcdf = "Gaussian", minSize = 10, maxSize = 500), verbose = FALSE)
paths <- c("HALLMARK_GLYCOLYSIS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_MTORC1_SIGNALING",
           "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_FATTY_ACID_METABOLISM",
           "HALLMARK_HYPOXIA", "HALLMARK_INFLAMMATORY_RESPONSE")
path_test <- function(sc, lab, cohort) {
  sid <- intersect(colnames(sc), lab$sample_id); cl <- lab$cluster[match(sid, lab$sample_id)]
  rbindlist(lapply(intersect(paths, rownames(sc)), function(pp) {
    a <- sc[pp, sid[cl == 1]]; b <- sc[pp, sid[cl == 2]]; w <- wilcox.test(a, b, exact = FALSE)
    data.table(cohort = cohort, pathway = pp, C1_mean = mean(a), C2_mean = mean(b),
               delta_C1_minus_C2 = mean(a) - mean(b), p_value = w$p.value)
  }))
}
pd <- rbindlist(list(path_test(g395sc, r395$labels, "GSE39582"), path_test(g175sc, r175$labels, "GSE17536")))
pd[, FDR := p.adjust(p_value, "BH"), by = cohort]
old_pd <- read_tsv(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage61/04_metabolic_heatmap/metabolic_pathway_direction_effects.tsv"))[cohort %chin% c("GSE39582", "GSE17536")]
cmp <- merge(old_pd[, .(cohort, pathway, old_delta = delta_C1_minus_C2, old_FDR = FDR_within_cohort)],
             pd[, .(cohort, pathway, new_delta = delta_C1_minus_C2, new_FDR = FDR)], by = c("cohort", "pathway"))
cmp[, direction_consistent := sign(old_delta) == sign(new_delta)]
fwrite(pd, file.path(out, "external_reprojection_Hallmark_differential.tsv"), sep = "\t")
fwrite(cmp, file.path(out, "old_vs_exclusion13_centroid_pathway_comparison.tsv"), sep = "\t")

# Fixed five-pathway equal-weight intensity under new labels (unadjusted label sensitivity).
five <- paths[1:5]
intensity <- function(sc, lab, cohort) {
  sid <- intersect(colnames(sc), lab$sample_id); cl <- lab$cluster[match(sid, lab$sample_id)]
  zz <- z_rows(sc[five, sid, drop = FALSE]); score <- colMeans(zz); w <- wilcox.test(score[cl == 1], score[cl == 2], exact = FALSE)
  data.table(cohort = cohort, n = length(sid), C1_mean = mean(score[cl == 1]), C2_mean = mean(score[cl == 2]),
             delta_C1_minus_C2 = mean(score[cl == 1]) - mean(score[cl == 2]), p_value = w$p.value)
}
inten <- rbindlist(list(intensity(g395sc, r395$labels, "GSE39582"), intensity(g175sc, r175$labels, "GSE17536")))
inten[, FDR := p.adjust(p_value, "BH")]
fwrite(inten, file.path(out, "external_reprojection_five_pathway_intensity.tsv"), sep = "\t")

draw_labels <- function() {
  op <- par(mfrow = c(1, 2), mar = c(5, 4, 3, 1)); on.exit(par(op))
  for (co in unique(samples$cohort)) {
    d <- samples[cohort == co]
    tab <- table(Old = d$old_cluster, New = d$new_cluster)
    barplot(tab, beside = FALSE, col = c("#6F9FB5", "#D28D68"), xlab = "New centroid label", ylab = "Samples", main = co)
    legend("topright", c("Old C1", "Old C2"), fill = c("#6F9FB5", "#D28D68"), bty = "n", cex = .8)
  }
}
save_plot(file.path(out, "old_vs_exclusion13_centroid_labels"), draw_labels, 9, 5)

audit <- data.table(item = c("TCGA_centroid_samples", "GSE39582_samples", "GSE17536_samples", "GSE39582_genes", "GSE17536_genes"),
                    value = c(length(ids), nrow(r395$sample), nrow(r175$sample), nrow(r395$z), nrow(r175$z)), status = "complete")
fwrite(audit, file.path(out, "reprojection_completion_audit.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
cat("External re-projection completed\n"); print(summary); print(margin_cmp); print(cmp); print(inten)
