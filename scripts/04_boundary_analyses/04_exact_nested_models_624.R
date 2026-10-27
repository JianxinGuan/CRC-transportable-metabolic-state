args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop('Usage: script patient.tsv marker.tsv.gz output_root')
patients <- read.delim(args[[1]], check.names=FALSE, stringsAsFactors=FALSE, comment.char='')
marker <- read.delim(gzfile(args[[2]]), check.names=FALSE, stringsAsFactors=FALSE, comment.char='')
out <- file.path(args[[3]], '04_CMS_MSI_rerun'); dir.create(out, recursive=TRUE, showWarnings=FALSE)
if(!'patient_barcode' %in% names(patients)) patients$patient_barcode <- patients$patient_barcode_fixed
marker <- marker[!duplicated(marker$sample_id),,drop=FALSE]
d <- merge(patients, marker, by='sample_id', all.x=FALSE, suffixes=c('_patient','_marker'))
if('cluster_patient' %in% names(d)) d$cluster <- d$cluster_patient
if('CMS_marker' %in% names(d)) d$CMS <- d$CMS_marker
if('project_marker' %in% names(d)) d$project <- d$project_marker
d$cluster <- factor(d$cluster); d$CMS <- factor(d$CMS, levels=paste0('CMS',1:4)); d$project <- factor(d$project)
write.table(d, file.path(out,'TCGA_624_exact_stage61_merged.tsv'), sep='\t', row.names=FALSE, quote=TRUE, na='')
nested_one <- function(y) {
 vars <- c(y,'CMS','acidosis_signed_score','proliferation_composite','microenv_proxy','project'); vars <- intersect(vars,names(d)); q <- d[complete.cases(d[,vars,drop=FALSE]),,drop=FALSE]
 base <- c('CMS','proliferation_composite','microenv_proxy','project'); base <- base[base %in% names(q)]; base <- base[vapply(q[,base,drop=FALSE], function(x) length(unique(x))>1, logical(1))]
 if(!length(base) || nrow(q)<100) return(data.frame(outcome=y,analysis='all_624_exact_stage61',n=nrow(q)))
 m0 <- lm(as.formula(paste(y,'~',paste(base,collapse=' + '))),q); m1 <- update(m0,.~.+acidosis_signed_score); aa <- anova(m0,m1); cf <- summary(m1)$coefficients['acidosis_signed_score',]; ci <- confint(m1,'acidosis_signed_score')
 data.frame(outcome=y,analysis='all_624_exact_stage61',n=nrow(q),state_beta=cf[1],SE=cf[2],CI_low=ci[1],CI_high=ci[2],state_p=cf[4],LRT_F=aa$F[2],LRT_p=aa$`Pr(>F)`[2],adjusted_R2_base=summary(m0)$adj.r.squared,adjusted_R2_full=summary(m1)$adj.r.squared,adjusted_R2_increment=summary(m1)$adj.r.squared-summary(m0)$adj.r.squared,partial_R2=(deviance(m0)-deviance(m1))/deviance(m0))
}
outmod <- do.call(rbind,lapply(c('acidosis_composite','metabolic_intensity'),nested_one)); outmod$FDR <- p.adjust(outmod$state_p,'BH'); write.table(outmod,file.path(out,'CMS_incremental_models_624_exact_stage61.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
within <- list(); ii <- 1L
for(y in c('acidosis_composite','metabolic_intensity')) for(cc in levels(d$CMS)) { ok <- as.character(d$CMS)==cc & complete.cases(d[,c(y,'cluster','acidosis_signed_score'),drop=FALSE]); q <- d[ok,,drop=FALSE]; if(nrow(q)>0 && length(unique(as.character(q$cluster)))==2) { p <- wilcox.test(q[[y]]~q$cluster,exact=FALSE)$p.value; c1 <- as.character(q$cluster)=='1'; c2 <- as.character(q$cluster)=='2'; within[[ii]] <- data.frame(outcome=y,CMS=cc,n=nrow(q),C1_n=sum(c1),C2_n=sum(c2),delta_C2_minus_C1=mean(q[c2,][[y]])-mean(q[c1,][[y]]),wilcox_p=p); ii<-ii+1L } }
within <- do.call(rbind,within); within$wilcox_FDR <- NA_real_; for(yy in unique(within$outcome)) within$wilcox_FDR[within$outcome==yy] <- p.adjust(within$wilcox_p[within$outcome==yy],'BH'); write.table(within,file.path(out,'within_CMS_effects_624_exact_stage61.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo_exact_stage61_624.txt'))
writeLines(c('Exact Stage 61 nested model replication on 624 deduplicated patients.','Formula: outcome ~ CMS + acidosis_signed_score + proliferation_composite + microenv_proxy + project.','Within-CMS effect uses mean C2-minus-C1 and Wilcoxon test, matching the original script.'),file.path(out,'exact_stage61_model_definition.txt'))
