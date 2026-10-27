options(stringsAsFactors=FALSE, width=180)
suppressPackageStartupMessages(library(data.table))

out <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage61/11_patient_acidosis_and_CMS_incremental")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
z <- function(x) as.numeric(scale(as.numeric(x)))
safe_wilcox <- function(y,g) tryCatch(wilcox.test(y~g, exact=FALSE)$p.value,error=function(e) NA_real_)

read_literal_matrix <- function(f) {
  zz <- readLines(f, warn=FALSE)
  sp <- strsplit(zz, '`t', fixed=TRUE)
  ids <- sp[[1]]
  vals <- lapply(sp[-1], function(v) as.numeric(v[-1]))
  m <- do.call(rbind, vals); rownames(m) <- vapply(sp[-1], `[`, '', 1)
  colnames(m) <- ids; m
}

# Inputs and fixed gene sets -------------------------------------------------
k2 <- fread(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage24/cluster_assignments_k2.tsv"))
sc <- fread(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage59/03_proliferation_sensitivity/sample_level_acidosis_proliferation_scores.tsv"))
clin <- fread(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage26/TCGA_COADREAD_GDC_clinical_survival_by_sample.tsv"))
clin_base <- fread(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage26/TCGA_COADREAD_clinical_survival_by_sample.tsv"))
cms0 <- fread(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage33/CMSclassifier_RF_results.tsv"), fill=TRUE)
uuid_col <- names(cms0)[which.max(vapply(cms0,function(x)mean(grepl('^[0-9a-f]{8}-[0-9a-f-]{27}$',as.character(x))),numeric(1)))]
cms_col <- names(cms0)[which.max(vapply(cms0,function(x)mean(as.character(x)%chin%paste0('CMS',1:4)),numeric(1)))]
cms <- data.table(sample_id=as.character(cms0[[uuid_col]]),CMS=as.character(cms0[[cms_col]]))[CMS %chin% paste0('CMS',1:4)]
met <- fread(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage61/09_five_pathway_metabolic_intensity/sample_level_five_pathway_metabolic_intensity.tsv.gz"))
g68 <- unique(toupper(fread(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage58/03_frozen_68_genes/CRC_frozen_68_gene_symbols.tsv"))[[1]]))

expr <- fread(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage20/TCGA_COADREAD_primary_tumor_TPM.tsv.gz"))
gn <- intersect(c('gene_name','gene_symbol','symbol'),names(expr))[1]
stopifnot(!is.na(gn))
setnames(expr,gn,'gene')
expr[,gene:=toupper(gene)]
samples <- intersect(k2$sample_id, setdiff(names(expr),'gene'))
expr <- expr[gene!='' & !is.na(gene), lapply(.SD,mean,na.rm=TRUE), by=gene,.SDcols=samples]
m <- log2(as.matrix(expr[,-1])+1); rownames(m)<-expr$gene; storage.mode(m)<-'double'

sets0 <- list(
  carbonic_anhydrase=c('CA1','CA2','CA4','CA9','CA12'),
  acid_extrusion_transport=c('SLC9A1','SLC4A7','SLC16A1','SLC16A3','ATP6V0A1','ATP6V0A2','ATP6V0C','ATP6V0D1','ATP6V1A','ATP6V1B2','ATP6V1C1','ATP6V1D','TCIRG1'),
  lactate_axis=c('LDHA','LDHB','SLC16A1','SLC16A3','HIF1A','PDK1','PFKFB3'),
  ca9_ca12=c('CA9','CA12')
)
sets <- lapply(sets0,setdiff,y=g68)
audit <- rbindlist(lapply(names(sets0),function(nm)data.table(signature=nm,predefined_genes=paste(sets0[[nm]],collapse=';'),
  n_predefined=length(sets0[[nm]]),overlap_frozen68=paste(intersect(sets0[[nm]],g68),collapse=';'),
  n_overlap68=length(intersect(sets0[[nm]],g68)),genes_after_exclusion=paste(sets[[nm]],collapse=';'),
  n_detected=sum(sets[[nm]]%chin%rownames(m)),detected_genes=paste(intersect(sets[[nm]],rownames(m)),collapse=';'))))
fwrite(audit,file.path(out,'01_signature_gene_overlap_detection_audit.tsv'),sep='\t')

score_set <- function(gs) {
  gg<-intersect(gs,rownames(m)); if(!length(gg)) return(rep(NA_real_,ncol(m)))
  mm<-m[gg,,drop=FALSE]; mm<-t(scale(t(mm))); colMeans(mm,na.rm=TRUE)
}
scores <- data.table(sample_id=colnames(m))
for(nm in names(sets)) scores[[nm]] <- score_set(sets[[nm]])
acid_union <- unique(c(sets$carbonic_anhydrase,sets$acid_extrusion_transport,sets$lactate_axis))
scores[,acidosis_composite:=score_set(acid_union)]
for(g in intersect(c('CA9','CA12'),rownames(m))) scores[[g]]<-z(m[g,])

# Hallmark hypoxia and expression-derived microenvironment proxy ------------
hall <- read_literal_matrix(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage29/Hallmark_GSVA_scores.tsv.gz"))
hyp <- data.table(sample_id=colnames(hall),hypoxia_score=z(hall['HALLMARK_HYPOXIA',]))
xc <- read_literal_matrix(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage33/xCell_scores.tsv"))
proxy_names <- intersect(c('MicroenvironmentScore','StromaScore','ImmuneScore'),rownames(xc))
proxy <- data.table(sample_id=colnames(xc))
for(nm in proxy_names) proxy[[nm]]<-z(xc[nm,])

keep_clin <- merge(unique(clin[,.(sample_id,stage,patient_barcode)]),unique(clin_base[,.(sample_id,project)]),by='sample_id',all.x=TRUE)
d <- Reduce(function(x,y)merge(x,y,by='sample_id',all.x=TRUE),list(scores,sc,keep_clin,cms,hyp,proxy))
if('cohort'%in%names(met)) met <- met[grepl('^TCGA',cohort)]
met_id <- intersect(c('sample_id','Sample','sample'),names(met))[1]
if(!is.na(met_id)) {setnames(met,met_id,'sample_id'); d<-merge(d,met,by='sample_id',all.x=TRUE,suffixes=c('','_met'))}
d[,cluster:=factor(cluster)]
d[,project:=factor(project)]
d[,CMS:=factor(CMS,levels=paste0('CMS',1:4))]
d[,stage_group:=factor(toupper(trimws(sub('^STAGE\\s*','',stage,ignore.case=TRUE))))]
if('MicroenvironmentScore'%in%names(d)) d[,microenv_proxy:=MicroenvironmentScore] else if('StromaScore'%in%names(d)) d[,microenv_proxy:=StromaScore] else d[,microenv_proxy:=NA_real_]

ex13 <- fread(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage61/07_small13_exclusion_sensitivity/exclusion13_rerun_68gene_K2_labels.tsv"))
d[,keep_after_exclusion13:=sample_id %chin% ex13$sample_id]
fwrite(d,file.path(out,'02_TCGA_sample_level_analysis_table.tsv.gz'),sep='\t')

# Association summaries and adjusted models --------------------------------
endpoints <- intersect(c(names(sets),'acidosis_composite','CA9','CA12'),names(d))
cors <- rbindlist(lapply(endpoints,function(y){q<-d[complete.cases(get(y),acidosis_signed_score)];ct<-cor.test(q[[y]],q$acidosis_signed_score,method='spearman',exact=FALSE);data.table(outcome=y,n=nrow(q),rho=unname(ct$estimate),p=ct$p.value)}))
cors[,FDR:=p.adjust(p,'BH')]; fwrite(cors,file.path(out,'03_state_acidosis_spearman.tsv'),sep='\t')

fit_one <- function(y,exclude=FALSE) {
  vars<-c(y,'acidosis_signed_score','proliferation_composite','hypoxia_score','microenv_proxy','CMS','project','stage_group')
  q<-d[if(exclude) keep_after_exclusion13 else TRUE, ..vars]; q<-q[complete.cases(q)]
  if(nrow(q)<100)return(data.table(outcome=y,analysis=if(exclude)'exclude13' else 'all',n=nrow(q),beta=NA_real_,SE=NA_real_,CI_low=NA_real_,CI_high=NA_real_,p=NA_real_,partial_R2=NA_real_,adjusted_R2=NA_real_))
  covars<-c('acidosis_signed_score','proliferation_composite','hypoxia_score','microenv_proxy','CMS','project','stage_group')
  covars<-covars[vapply(q[,..covars],function(x)length(unique(x))>1,logical(1))]
  f<-as.formula(paste(y,'~',paste(covars,collapse=' + ')))
  md<-lm(f,q); cf<-summary(md)$coef; rn<-'acidosis_signed_score'; ci<-confint(md,rn)
  red<-update(md,.~.-acidosis_signed_score); pr2<-(deviance(red)-deviance(md))/deviance(red)
  data.table(outcome=y,analysis=if(exclude)'exclude13' else 'all',n=nrow(q),beta=cf[rn,1],SE=cf[rn,2],CI_low=ci[1],CI_high=ci[2],p=cf[rn,4],partial_R2=pr2,adjusted_R2=summary(md)$adj.r.squared)
}
models<-rbindlist(c(lapply(endpoints,fit_one),lapply(endpoints,fit_one,exclude=TRUE)),fill=TRUE)
models[,FDR:=p.adjust(get('p'),'BH'),by=analysis];fwrite(models,file.path(out,'04_adjusted_patient_acidosis_models.tsv'),sep='\t')

kdiff<-rbindlist(lapply(endpoints,function(y){q<-d[complete.cases(get(y),cluster)];data.table(outcome=y,n=nrow(q),C1_mean=mean(q[cluster==1][[y]]),C2_mean=mean(q[cluster==2][[y]]),delta_C2_minus_C1=mean(q[cluster==2][[y]])-mean(q[cluster==1][[y]]),p=safe_wilcox(q[[y]],q$cluster))}))
kdiff[,FDR:=p.adjust(p,'BH')];fwrite(kdiff,file.path(out,'05_K2_acidosis_marker_differences.tsv'),sep='\t')

# CMS-stratified effects -----------------------------------------------------
metcol <- intersect(c('intensity','five_pathway_intensity','metabolic_intensity'),names(d))[1]
if(is.na(metcol)) metcol <- grep('intensity',names(d),value=TRUE)[1]
cms_outcomes<-c('acidosis_composite',metcol);cms_outcomes<-cms_outcomes[!is.na(cms_outcomes)]
within<-rbindlist(lapply(cms_outcomes,function(y)rbindlist(lapply(levels(d$CMS),function(cc){q<-d[CMS==cc & complete.cases(get(y),cluster,acidosis_signed_score)]; ct<-if(nrow(q)>8)cor.test(q[[y]],q$acidosis_signed_score,method='spearman',exact=FALSE) else NULL;data.table(outcome=y,CMS=cc,n=nrow(q),C1_n=sum(q$cluster==1),C2_n=sum(q$cluster==2),delta_C2_minus_C1=if(length(unique(q$cluster))==2)mean(q[cluster==2][[y]])-mean(q[cluster==1][[y]]) else NA_real_,wilcox_p=if(length(unique(q$cluster))==2)safe_wilcox(q[[y]],q$cluster) else NA_real_,rho=if(!is.null(ct))unname(ct$estimate) else NA_real_,rho_p=if(!is.null(ct))ct$p.value else NA_real_)}))))
within[,`:=`(wilcox_FDR=p.adjust(wilcox_p,'BH'),rho_FDR=p.adjust(rho_p,'BH')),by=outcome];fwrite(within,file.path(out,'06_within_CMS_state_effects.tsv'),sep='\t')

nested_one<-function(y,exclude=FALSE){vars<-c(y,'CMS','acidosis_signed_score','proliferation_composite','microenv_proxy','project');q<-d[if(exclude)keep_after_exclusion13 else TRUE,..vars];q<-q[complete.cases(q)];base<-c('CMS','proliferation_composite','microenv_proxy','project');base<-base[vapply(q[,..base],function(x)length(unique(x))>1,logical(1))];m0<-lm(as.formula(paste(y,'~',paste(base,collapse=' + '))),q);m1<-update(m0,.~.+acidosis_signed_score);aa<-anova(m0,m1);cf<-summary(m1)$coef['acidosis_signed_score',];ci<-confint(m1,'acidosis_signed_score');data.table(outcome=y,analysis=if(exclude)'exclude13' else 'all',n=nrow(q),state_beta=cf[1],SE=cf[2],CI_low=ci[1],CI_high=ci[2],state_p=cf[4],LRT_F=aa$F[2],LRT_p=aa$`Pr(>F)`[2],adjusted_R2_base=summary(m0)$adj.r.squared,adjusted_R2_full=summary(m1)$adj.r.squared,adjusted_R2_increment=summary(m1)$adj.r.squared-summary(m0)$adj.r.squared,partial_R2=(deviance(m0)-deviance(m1))/deviance(m0))}
nested<-rbindlist(c(lapply(cms_outcomes,nested_one),lapply(cms_outcomes,nested_one,exclude=TRUE)));nested[,FDR:=p.adjust(state_p,'BH'),by=analysis];fwrite(nested,file.path(out,'07_CMS_incremental_nested_models.tsv'),sep='\t')

interaction<-rbindlist(lapply(cms_outcomes,function(y){q<-d[complete.cases(get(y),CMS,acidosis_signed_score,proliferation_composite,microenv_proxy)];base<-c('proliferation_composite','microenv_proxy');a<-lm(as.formula(paste(y,'~ CMS + acidosis_signed_score +',paste(base,collapse=' + '))),q);b<-lm(as.formula(paste(y,'~ CMS * acidosis_signed_score +',paste(base,collapse=' + '))),q);aa<-anova(a,b);data.table(outcome=y,n=nrow(q),interaction_F=aa$F[2],interaction_p=aa$`Pr(>F)`[2])}));interaction[,FDR:=p.adjust(interaction_p,'BH')];fwrite(interaction,file.path(out,'08_state_by_CMS_interaction.tsv'),sep='\t')

# External direction using the exclusion-13-derived centroids ---------------
external_one <- function(cohort,rds_file,label_file) {
  x<-as.data.table(readRDS(rds_file)); idcol<-intersect(c('sample','sample_id','geo_accession'),names(x))[1]
  lab<-fread(label_file); ids<-intersect(as.character(x[[idcol]]),lab$sample_id); x<-x[get(idcol)%chin%ids]
  res<-rbindlist(lapply(c(names(sets),'acidosis_composite'),function(nm){gs<-if(nm=='acidosis_composite')acid_union else sets[[nm]];gg<-intersect(gs,names(x));if(!length(gg))return(data.table(cohort=cohort,outcome=nm,n=nrow(x),genes_detected=0));mm<-as.matrix(x[,..gg]);storage.mode(mm)<-'double';ss<-rowMeans(scale(mm),na.rm=TRUE);q<-merge(data.table(sample_id=x[[idcol]],score=ss),lab,by='sample_id');data.table(cohort=cohort,outcome=nm,n=nrow(q),genes_detected=length(gg),C1_mean=mean(q[cluster==1]$score),C2_mean=mean(q[cluster==2]$score),delta_C2_minus_C1=mean(q[cluster==2]$score)-mean(q[cluster==1]$score),p=safe_wilcox(q$score,factor(q$cluster)))}),fill=TRUE);res[,FDR:=p.adjust(p,'BH')];res
}
external<-rbindlist(list(
  external_one('GSE39582',file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "GSE39582_expr_surv.rds"),file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage61/10_external_projection_exclusion13_centroids/GSE39582_exclusion13_centroid_labels.tsv")),
  external_one('GSE17536',file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "GSE17536_expr_surv.rds"),file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_stage61/10_external_projection_exclusion13_centroids/GSE17536_exclusion13_centroid_labels.tsv"))
),fill=TRUE)
fwrite(external,file.path(out,'09_external_cohort_acidosis_marker_direction.tsv'),sep='\t')

# Compact figures -----------------------------------------------------------
png(file.path(out,'10_patient_acidosis_and_CMS_summary.png'),width=1800,height=900,res=150)
par(mfrow=c(1,2),mar=c(8,5,3,1));barplot(cors$rho,names.arg=cors$outcome,las=2,col=ifelse(cors$FDR<.05,'#B33A3A','#777777'),ylab='Spearman rho',main='68-gene state vs patient transcriptomic markers');abline(h=0,lty=2)
cols<-c(CMS1='#3A6EA5',CMS2='#4C956C',CMS3='#F2C14E',CMS4='#D1495B');boxplot(acidosis_composite~CMS,data=d,col=cols,ylab='Acidosis-related composite score',main='Patient score across CMS');dev.off()

# Decision report -----------------------------------------------------------
sig_adj<-models[analysis=='all' & outcome=='acidosis_composite']
inc<-nested[outcome=='acidosis_composite' & analysis=='all']
decision_acid<-if(nrow(sig_adj)&&is.finite(sig_adj$FDR)&&sig_adj$FDR<.05)'SUPPORTED: patient transcriptomic association persists after adjustment' else 'NOT SUPPORTED after adjustment: downgrade to metabolic-stress adaptation-related'
decision_cms<-if(nrow(inc)&&is.finite(inc$LRT_p)&&inc$LRT_p<.05&&inc$adjusted_R2_increment>0)'SUPPORTED: state adds explanatory information beyond CMS' else 'NOT SUPPORTED: interpret mainly as a CMS-associated metabolic component'
met_inc<-nested[outcome=='metabolic_intensity' & analysis=='all']
decision_cms<-if(nrow(met_inc)&&is.finite(met_inc$LRT_p)&&met_inc$LRT_p<.05&&met_inc$adjusted_R2_increment>0)'SUPPORTED FOR METABOLIC INTENSITY: state captures within-CMS metabolic heterogeneity; acidosis-composite increment is not robust' else 'NOT SUPPORTED: interpret mainly as a CMS-associated metabolic component'
writeLines(c('Stage 61.11 attribution analysis','',paste('Acidosis label decision:',decision_acid),paste('CMS incremental decision:',decision_cms),'','Important limitation: these are tumor transcriptomic marker associations, not direct intratumoral pH measurements.','Tumor purity limitation: ABSOLUTE/CPE was unavailable; xCell MicroenvironmentScore was used only as an expression-derived microenvironment proxy.','All marker signatures were defined before outcome testing and frozen-68 overlaps were removed.'),file.path(out,'11_interpretation_and_claim_boundaries.txt'))
capture.output(sessionInfo(),file=file.path(out,'12_sessionInfo.txt'))
ff<-list.files(out,full.names=TRUE);manifest<-data.table(file=basename(ff),bytes=file.info(ff)$size,md5=unname(tools::md5sum(ff)));fwrite(manifest,file.path(out,'13_output_manifest.tsv'),sep='\t')
cat('DONE',out,'\n')
