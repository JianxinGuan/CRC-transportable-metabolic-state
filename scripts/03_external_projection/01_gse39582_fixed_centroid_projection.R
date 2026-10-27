suppressPackageStartupMessages({
  library(data.table)
  library(cluster)
  library(survival)
  library(GSVA)
})

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name) { hit <- match(name, args); if (is.na(hit) || hit == length(args)) stop("Missing ", name); args[[hit + 1L]] }
out <- arg("--out")
tcga_z_file <- arg("--tcga-zscore")
tcga_label_file <- arg("--tcga-labels")
exclude13_label_file <- arg("--exclude13-labels")
g395_rds_file <- arg("--gse39582-rds")
g175_rds_file <- arg("--gse17536-rds")
g395_original_label_file <- arg("--gse39582-labels")
g175_original_label_file <- arg("--gse17536-labels")
g395_score_file <- arg("--gse39582-hallmark-scores")
tcga_score_file <- arg("--tcga-hallmark-scores")
gmt_file <- arg("--hallmark-gmt")
proliferation_file <- arg("--tcga-proliferation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dirs <- file.path(out, c("original_projection", "exclusion13_projection", "pathways", "intensity", "balance", "acidification"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

z_rows <- function(matrix) { result <- t(scale(t(matrix))); result[!is.finite(result)] <- 0; result }
ari <- function(x, y) {
  tab <- table(x, y); n <- sum(tab); c2 <- function(z) z * (z - 1) / 2
  a <- sum(c2(tab)); b <- sum(c2(rowSums(tab))); c <- sum(c2(colSums(tab)))
  expected <- b * c / c2(n); denominator <- (b + c) / 2 - expected
  if (denominator == 0) NA_real_ else (a - expected) / denominator
}
nmi <- function(x, y) {
  tab <- table(x, y); p <- tab / sum(tab); px <- rowSums(p); py <- colSums(p); expected <- outer(px, py); nonzero <- p > 0
  mutual <- sum(p[nonzero] * log(p[nonzero] / expected[nonzero]))
  mutual / sqrt((-sum(px[px > 0] * log(px[px > 0]))) * (-sum(py[py > 0] * log(py[py > 0]))))
}

tcga <- fread(tcga_z_file)
genes <- sub("^\\ufeff", "", tcga[[1L]])
tcga_matrix <- as.matrix(tcga[, -1L]); rownames(tcga_matrix) <- genes; storage.mode(tcga_matrix) <- "double"
tcga_labels <- fread(tcga_label_file)
tcga_cluster <- tcga_labels$cluster[match(colnames(tcga_matrix), tcga_labels$sample_id)]
stopifnot(!anyNA(tcga_cluster))
original_centroids <- sapply(1:2, function(group) rowMeans(tcga_matrix[, tcga_cluster == group, drop = FALSE]))
exclude_labels <- fread(exclude13_label_file)
exclude_ids <- intersect(exclude_labels$sample_id, colnames(tcga_matrix))
exclude_labels <- exclude_labels[match(exclude_ids, sample_id)]
exclude_matrix <- tcga_matrix[, exclude_ids, drop = FALSE]
exclude_centroids <- sapply(1:2, function(group) rowMeans(exclude_matrix[, exclude_labels$exclusion13_rerun_K2 == group, drop = FALSE]))

g395 <- as.data.table(readRDS(g395_rds_file))
g175 <- as.data.table(readRDS(g175_rds_file))
make_z <- function(x, requested) {
  available <- intersect(requested, names(x))
  matrix <- t(as.matrix(x[, ..available])); colnames(matrix) <- as.character(x$sample)
  list(z = z_rows(matrix), genes = available)
}
z395 <- make_z(g395, genes)
z175 <- make_z(g175, genes)
stopifnot(length(z395$genes) == 68L, length(z175$genes) == 68L)

projection_metrics <- function(z, labels, centroids, cohort) {
  ids <- intersect(colnames(z), labels$sample_id); z <- z[, ids, drop = FALSE]
  assigned <- labels$cluster[match(ids, labels$sample_id)]
  silhouette_width <- silhouette(assigned, as.dist(1 - cor(z)))[, "sil_width"]
  cor1 <- apply(z, 2, cor, centroids[rownames(z), 1L]); cor2 <- apply(z, 2, cor, centroids[rownames(z), 2L])
  external_centroids <- sapply(1:2, function(group) rowMeans(z[, assigned == group, drop = FALSE]))
  sample <- data.table(cohort = cohort, sample_id = ids, cluster = assigned, silhouette = silhouette_width,
                       cor_TCGA_C1 = cor1, cor_TCGA_C2 = cor2, assignment_margin = abs(cor1 - cor2))
  summary <- data.table(cohort = cohort, n = length(ids), genes = nrow(z), cluster1_n = sum(assigned == 1L),
                        cluster2_n = sum(assigned == 2L), cluster1_percent = 100 * mean(assigned == 1L),
                        cluster2_percent = 100 * mean(assigned == 2L),
                        imbalance_ratio_larger_to_smaller = max(table(assigned)) / min(table(assigned)),
                        silhouette_mean = mean(silhouette_width), silhouette_median = median(silhouette_width),
                        negative_silhouette_percent = 100 * mean(silhouette_width < 0),
                        assignment_margin_median = median(abs(cor1 - cor2)),
                        centroid_distance_C1 = 1 - cor(external_centroids[, 1L], centroids[rownames(z), 1L]),
                        centroid_distance_C2 = 1 - cor(external_centroids[, 2L], centroids[rownames(z), 2L]))
  list(sample = sample, summary = summary)
}
labels395 <- fread(g395_original_label_file)
labels175 <- fread(g175_original_label_file)
original395 <- projection_metrics(z395$z, labels395, original_centroids, "GSE39582")
original175 <- projection_metrics(z175$z, labels175, original_centroids, "GSE17536")
original_samples <- rbindlist(list(original395$sample, original175$sample))
original_summary <- rbindlist(list(original395$summary, original175$summary))
fwrite(original_samples, file.path(dirs[1L], "sample_level_projection_quality.tsv.gz"), sep = "\t")
fwrite(original_summary, file.path(dirs[1L], "projection_quality_and_subtype_sizes.tsv"), sep = "\t")

project_new <- function(x, z, old, cohort) {
  cor1 <- apply(z, 2, cor, exclude_centroids[rownames(z), 1L]); cor2 <- apply(z, 2, cor, exclude_centroids[rownames(z), 2L])
  prediction <- ifelse(cor2 > cor1, 2L, 1L); ids <- colnames(z)
  old_cluster <- old$cluster[match(ids, old$sample_id)]; old_margin <- old$margin[match(ids, old$sample_id)]
  sil <- silhouette(prediction, as.dist(1 - cor(z)))[, "sil_width"]
  external_centroids <- sapply(1:2, function(group) rowMeans(z[, prediction == group, drop = FALSE]))
  os_status <- as.numeric(as.character(x$OS.status)); os_time <- as.numeric(x$OS.time)
  complete <- is.finite(os_status) & is.finite(os_time) & os_time > 0
  surv <- survdiff(Surv(os_time[complete], os_status[complete]) ~ prediction[complete])
  sample <- data.table(cohort = cohort, sample_id = ids, old_cluster = old_cluster, new_cluster = prediction,
                       changed = prediction != old_cluster, old_margin = old_margin, new_cor_C1 = cor1,
                       new_cor_C2 = cor2, new_margin = abs(cor1 - cor2), silhouette = sil,
                       OS_status = os_status, OS_time = os_time)
  summary <- data.table(cohort = cohort, n = length(ids), genes = nrow(z), old_cluster1_n = sum(old_cluster == 1L),
                        old_cluster2_n = sum(old_cluster == 2L), new_cluster1_n = sum(prediction == 1L),
                        new_cluster2_n = sum(prediction == 2L), changed_n = sum(prediction != old_cluster),
                        agreement = mean(prediction == old_cluster), ARI = ari(prediction, old_cluster), NMI = nmi(prediction, old_cluster),
                        new_silhouette_mean = mean(sil), new_silhouette_median = median(sil),
                        new_negative_silhouette_percent = 100 * mean(sil < 0), new_margin_median = median(abs(cor1-cor2)),
                        centroid_distance_C1 = 1-cor(external_centroids[,1L],exclude_centroids[rownames(z),1L]),
                        centroid_distance_C2 = 1-cor(external_centroids[,2L],exclude_centroids[rownames(z),2L]),
                        new_OS_logrank_p = pchisq(surv$chisq, 1, lower.tail = FALSE))
  list(sample = sample, summary = summary, labels = data.table(sample_id = ids, cluster = prediction))
}
new395 <- project_new(g395, z395$z, labels395, "GSE39582")
new175 <- project_new(g175, z175$z, labels175, "GSE17536")
new_samples <- rbindlist(list(new395$sample, new175$sample)); new_summary <- rbindlist(list(new395$summary, new175$summary))
fwrite(new_samples, file.path(dirs[2L], "external_reprojection_sample_level.tsv.gz"), sep = "\t")
fwrite(new_summary, file.path(dirs[2L], "external_reprojection_quality_summary.tsv"), sep = "\t")
fwrite(new395$labels, file.path(dirs[2L], "GSE39582_exclusion13_centroid_labels.tsv"), sep = "\t")
fwrite(new175$labels, file.path(dirs[2L], "GSE17536_exclusion13_centroid_labels.tsv"), sep = "\t")
margin <- rbindlist(lapply(unique(new_samples$cohort), function(cohort_value) {
  d <- new_samples[new_samples$cohort == cohort_value]; test <- wilcox.test(d$old_margin[d$changed], d$old_margin[!d$changed], exact = FALSE)
  data.table(cohort = cohort_value, changed_n = sum(d$changed), unchanged_n = sum(!d$changed),
             changed_old_margin_median = median(d$old_margin[d$changed]),
             unchanged_old_margin_median = median(d$old_margin[!d$changed]), p_value = test$p.value)
}))
margin[, FDR := p.adjust(p_value, "BH")]
fwrite(margin, file.path(dirs[2L], "changed_sample_margin_comparison.tsv"), sep = "\t")

read_score <- function(path) {
  d <- fread(path)
  if (ncol(d) == 1L) {
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path) else file(path)
    x <- readLines(con, warn = FALSE)
    close(con)
    x <- gsub("`t", "\t", x, fixed = TRUE)
    d <- fread(text = paste(x, collapse = "\n"), sep = "\t")
  }
  m <- as.matrix(d[, -1L])
  rownames(m) <- d[[1L]]
  storage.mode(m) <- "double"
  m
}
score395 <- read_score(g395_score_file)
score_tcga <- read_score(tcga_score_file)
gmt <- strsplit(readLines(gmt_file, warn = FALSE), "\t")
sets <- setNames(lapply(gmt, function(value) value[-c(1L,2L)]), vapply(gmt, `[[`, character(1), 1L))
g175_genes <- setdiff(names(g175), c("sample","OS.status","OS.time")); matrix175 <- t(as.matrix(g175[, ..g175_genes])); colnames(matrix175) <- as.character(g175$sample)
sets175 <- lapply(sets, intersect, rownames(matrix175)); sets175 <- sets175[lengths(sets175) >= 10L]
score175 <- gsva(gsvaParam(matrix175, sets175, kcdf="Gaussian", minSize=10, maxSize=500), verbose=FALSE)
paths <- c("HALLMARK_GLYCOLYSIS","HALLMARK_MTORC1_SIGNALING","HALLMARK_OXIDATIVE_PHOSPHORYLATION",
           "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY","HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_HYPOXIA","HALLMARK_INFLAMMATORY_RESPONSE")
diff_table <- function(scores, labels, cohort) {
  ids <- intersect(colnames(scores), labels$sample_id); groups <- labels$cluster[match(ids, labels$sample_id)]
  result <- rbindlist(lapply(intersect(paths, rownames(scores)), function(pathway) {
    a <- scores[pathway, ids[groups==1L]]; b <- scores[pathway, ids[groups==2L]]; test <- wilcox.test(a,b,exact=FALSE)
    data.table(cohort=cohort,pathway=pathway,C1_mean=mean(a),C2_mean=mean(b),delta_C1_minus_C2=mean(a)-mean(b),p_value=test$p.value)
  })); result[,FDR_within_cohort:=p.adjust(p_value,"BH")]; result[,display_scaled_delta:=delta_C1_minus_C2/max(abs(delta_C1_minus_C2))]; result
}
pathway_results <- rbindlist(list(diff_table(score_tcga,tcga_labels,"TCGA_COADREAD"),diff_table(score395,labels395,"GSE39582"),diff_table(score175,labels175,"GSE17536")))
fwrite(pathway_results,file.path(dirs[3L],"metabolic_pathway_direction_effects.tsv"),sep="\t")

e2g2 <- unique(c(sets[["HALLMARK_E2F_TARGETS"]],sets[["HALLMARK_G2M_CHECKPOINT"]]))
prolif_external <- function(x) { ids<-as.character(x$sample); gg<-intersect(e2g2,names(x)); m<-t(as.matrix(x[,..gg])); colnames(m)<-ids; colMeans(z_rows(m)) }
tcga_prolif_dt <- fread(proliferation_file); tcga_prolif <- setNames(tcga_prolif_dt$proliferation_composite,tcga_prolif_dt$sample_id)
prolif395 <- prolif_external(g395); prolif175 <- prolif_external(g175)
five <- paths[1:5]
intensity_analysis <- function(scores, labels, proliferation, cohort) {
  ids <- Reduce(intersect,list(colnames(scores),labels$sample_id,names(proliferation))); standardized <- z_rows(scores[five,ids,drop=FALSE]); value<-colMeans(standardized)
  d<-data.table(cohort=cohort,sample_id=ids,cluster=factor(labels$cluster[match(ids,labels$sample_id)]),metabolic_intensity=value,proliferation=proliferation[ids])
  fit0<-lm(metabolic_intensity~cluster,d);fit1<-lm(metabolic_intensity~cluster+proliferation,d);c0<-summary(fit0)$coef["cluster2",];c1<-summary(fit1)$coef["cluster2",]
  result<-data.table(cohort=cohort,n=nrow(d),C1_mean=mean(d[cluster==1,metabolic_intensity]),C2_mean=mean(d[cluster==2,metabolic_intensity]),
                     delta_C1_minus_C2=mean(d[cluster==1,metabolic_intensity])-mean(d[cluster==2,metabolic_intensity]),
                     unadjusted_beta_C2_vs_C1=c0["Estimate"],unadjusted_p=c0["Pr(>|t|)"],adjusted_beta_C2_vs_C1=c1["Estimate"],adjusted_p=c1["Pr(>|t|)"],
                     proliferation_beta=summary(fit1)$coef["proliferation","Estimate"])
  list(sample=d,result=result)
}
i_tcga<-intensity_analysis(score_tcga,tcga_labels,tcga_prolif,"TCGA_COADREAD");i395<-intensity_analysis(score395,labels395,prolif395,"GSE39582");i175<-intensity_analysis(score175,labels175,prolif175,"GSE17536")
intensity_results<-rbindlist(list(i_tcga$result,i395$result,i175$result));intensity_results[,`:=`(unadjusted_FDR=p.adjust(unadjusted_p,"BH"),adjusted_FDR=p.adjust(adjusted_p,"BH"))]
intensity_samples<-rbindlist(list(i_tcga$sample,i395$sample,i175$sample));fwrite(intensity_results,file.path(dirs[4L],"five_pathway_metabolic_intensity_K2_results.tsv"),sep="\t");fwrite(intensity_samples,file.path(dirs[4L],"sample_level_five_pathway_metabolic_intensity.tsv.gz"),sep="\t")

balance_one <- function(scores,labels,proliferation,cohort) {
  ids<-Reduce(intersect,list(colnames(scores),labels$sample_id,names(proliferation)));d<-data.table(sample_id=ids,cluster=factor(labels$cluster[match(ids,labels$sample_id)]),glycolysis=scores["HALLMARK_GLYCOLYSIS",ids],OXPHOS=scores["HALLMARK_OXIDATIVE_PHOSPHORYLATION",ids],proliferation=proliferation[ids]);d[,balance_GLY_minus_OXPHOS:=as.numeric(scale(glycolysis))-as.numeric(scale(OXPHOS))];d[,coordinated_metabolic_activity:=(as.numeric(scale(glycolysis))+as.numeric(scale(OXPHOS)))/2]
  get_result<-function(y,outcome){f0<-lm(y~cluster,d);f1<-lm(y~cluster+proliferation,d);c0<-summary(f0)$coef["cluster2",];c1<-summary(f1)$coef["cluster2",];data.table(cohort=cohort,outcome=outcome,n=nrow(d),C1_median=median(y[d$cluster==1]),C2_median=median(y[d$cluster==2]),delta_C1_minus_C2=median(y[d$cluster==1])-median(y[d$cluster==2]),unadjusted_beta_C2_vs_C1=c0["Estimate"],unadjusted_p=c0["Pr(>|t|)"],adjusted_beta_C2_vs_C1=c1["Estimate"],adjusted_p=c1["Pr(>|t|)"],proliferation_beta=summary(f1)$coef["proliferation","Estimate"])}
  rbind(get_result(d$balance_GLY_minus_OXPHOS,"GLY_minus_OXPHOS_balance"),get_result(d$coordinated_metabolic_activity,"GLY_OXPHOS_coordinated_activity"))
}
balance_results<-rbindlist(list(balance_one(score_tcga,tcga_labels,tcga_prolif,"TCGA_COADREAD"),balance_one(score395,labels395,prolif395,"GSE39582"),balance_one(score175,labels175,prolif175,"GSE17536")));balance_results[,adjusted_FDR:=p.adjust(adjusted_p,"BH"),by=outcome];fwrite(balance_results,file.path(dirs[5L],"glycolysis_OXPHOS_balance_sensitivity.tsv"),sep="\t")

acid_sets<-list(carbonic_anhydrase=c("CA1","CA2","CA4","CA9","CA12"),acid_extrusion_transport=c("SLC9A1","SLC4A7","SLC16A1","SLC16A3","ATP6V0A1","ATP6V0A2","ATP6V0C","ATP6V0D1","ATP6V1A","ATP6V1B2","ATP6V1C1","ATP6V1D","TCIRG1"),lactate_axis=c("LDHA","LDHB","SLC16A1","SLC16A3","HIF1A","PDK1","PFKFB3"),ca9_ca12=c("CA9","CA12"));acid_union<-unique(c(acid_sets$carbonic_anhydrase,acid_sets$acid_extrusion_transport,acid_sets$lactate_axis))
safe_wilcox<-function(value,group) wilcox.test(value~group,exact=FALSE)$p.value
external_acid<-function(cohort,x,labels){result<-rbindlist(lapply(c(names(acid_sets),"acidosis_composite"),function(name){set<-if(name=="acidosis_composite")acid_union else acid_sets[[name]];detected<-intersect(set,names(x));matrix<-as.matrix(x[,..detected]);storage.mode(matrix)<-"double";score<-rowMeans(scale(matrix),na.rm=TRUE);q<-merge(data.table(sample_id=as.character(x$sample),score=score),labels,by="sample_id");data.table(cohort=cohort,outcome=name,n=nrow(q),genes_detected=length(detected),C1_mean=mean(q[cluster==1]$score),C2_mean=mean(q[cluster==2]$score),delta_C2_minus_C1=mean(q[cluster==2]$score)-mean(q[cluster==1]$score),p=safe_wilcox(q$score,factor(q$cluster)))}));result[,FDR:=p.adjust(p,"BH")];result}
acid_results<-rbindlist(list(external_acid("GSE39582",g395,new395$labels),external_acid("GSE17536",g175,new175$labels)));fwrite(acid_results,file.path(dirs[6L],"external_cohort_acidosis_marker_direction.tsv"),sep="\t")
writeLines(capture.output(sessionInfo()),file.path(out,"R_sessionInfo.txt"));cat("Corrected affected stage61 modules completed\n");print(original_summary);print(new_summary);print(margin);print(intensity_results);print(acid_results)
