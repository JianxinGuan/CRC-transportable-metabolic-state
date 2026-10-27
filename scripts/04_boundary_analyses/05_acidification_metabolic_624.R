args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=5L) stop('Usage: script patient.tsv marker.tsv.gz hallmark.txt intensity.tsv.gz output_root')
patient_file<-args[1]; marker_file<-args[2]; hallmark_file<-args[3]; intensity_file<-args[4]; root<-args[5]
out<-file.path(root,'05_acidification_metabolic_rerun');dir.create(out,recursive=TRUE,showWarnings=FALSE)
read_literal_matrix<-function(f){zz<-readLines(f,warn=FALSE);sp<-strsplit(zz,'`t',fixed=TRUE);ids<-sp[[1]];vals<-lapply(sp[-1],function(v)as.numeric(v[-1]));m<-do.call(rbind,vals);rownames(m)<-vapply(sp[-1],`[`,'',1);colnames(m)<-ids;m}
patients<-read.delim(patient_file,check.names=FALSE,stringsAsFactors=FALSE,comment.char='')
marker<-read.delim(gzfile(marker_file),check.names=FALSE,stringsAsFactors=FALSE,comment.char='')
inten<-read.delim(gzfile(intensity_file),check.names=FALSE,stringsAsFactors=FALSE,comment.char='')
hall<-read_literal_matrix(hallmark_file)
marker<-marker[!duplicated(marker$sample_id),]; tc_int<-inten[grepl('^TCGA',inten$cohort),];tc_int<-tc_int[!duplicated(tc_int$sample_id),]
d<-merge(patients,marker,by='sample_id',all.x=FALSE,suffixes=c('_patient','_marker'));d<-merge(d,tc_int[,c('sample_id','metabolic_intensity','proliferation')],by='sample_id',all.x=TRUE,suffixes=c('','_intensity'))
d$cluster<-factor(d$cluster_patient,levels=c(1,2));d$CMS<-factor(d$CMS,levels=paste0('CMS',1:4));d$project<-factor(d$project);d$stage_group<-factor(d$stage_group)
five<-c('HALLMARK_GLYCOLYSIS','HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_MTORC1_SIGNALING','HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY','HALLMARK_FATTY_ACID_METABOLISM')
for(nm in five)d[[nm]]<-as.numeric(hall[nm,d$sample_id])
write.table(d,file.path(out,'TCGA_624_acidification_metabolic_analysis_table.tsv'),sep='\t',row.names=FALSE,quote=TRUE,na='')

## Categorical C2-C1 differences for acidification composite and lactate axis.
catres<-do.call(rbind,lapply(c('acidosis_composite','lactate_axis'),function(y){q<-d[complete.cases(d[,c(y,'cluster')]),];w<-wilcox.test(q[[y]]~q$cluster,exact=FALSE);data.frame(outcome=y,n=nrow(q),C1_n=sum(q$cluster==1),C2_n=sum(q$cluster==2),C1_mean=mean(q[q$cluster==1,y]),C2_mean=mean(q[q$cluster==2,y]),delta_C2_minus_C1=mean(q[q$cluster==2,y])-mean(q[q$cluster==1,y]),wilcox_p=w$p.value)}));catres$FDR<-p.adjust(catres$wilcox_p,'BH');write.table(catres,file.path(out,'K2_acidification_lactate_differences_624.tsv'),sep='\t',row.names=FALSE,quote=TRUE)

## Exact original fully adjusted continuous-state models; FDR across the original seven outcomes.
endpoints<-c('carbonic_anhydrase','acid_extrusion_transport','lactate_axis','ca9_ca12','acidosis_composite','CA9','CA12')
fit_one<-function(y){vars<-c(y,'acidosis_signed_score','proliferation_composite','hypoxia_score','microenv_proxy','CMS','project','stage_group');q<-d[complete.cases(d[,vars]),vars,drop=FALSE];covars<-vars[-1];covars<-covars[vapply(q[,covars,drop=FALSE],function(x)length(unique(x))>1,logical(1))];md<-lm(as.formula(paste(y,'~',paste(covars,collapse=' + '))),q);cf<-summary(md)$coef['acidosis_signed_score',];ci<-confint(md,'acidosis_signed_score');red<-update(md,.~.-acidosis_signed_score);data.frame(outcome=y,n=nrow(q),beta=cf[1],SE=cf[2],CI_low=ci[1],CI_high=ci[2],p=cf[4],partial_R2=(deviance(red)-deviance(md))/deviance(red),adjusted_R2=summary(md)$adj.r.squared)}
adj<-do.call(rbind,lapply(endpoints,fit_one));adj$FDR<-p.adjust(adj$p,'BH');write.table(adj,file.path(out,'adjusted_acidification_models_624_exact_stage61.tsv'),sep='\t',row.names=FALSE,quote=TRUE)

## Five individual Hallmark pathway effects, BH corrected within the five-pathway family.
pathres<-do.call(rbind,lapply(five,function(y){q<-d[complete.cases(d[,c(y,'cluster')]),];fit<-lm(q[[y]]~q$cluster);cf<-summary(fit)$coef[2,];data.frame(pathway=y,n=nrow(q),C1_mean=mean(q[q$cluster==1,y]),C2_mean=mean(q[q$cluster==2,y]),beta_C2_vs_C1=cf[1],p=cf[4])}));pathres$FDR<-p.adjust(pathres$p,'BH');write.table(pathres,file.path(out,'five_pathway_effects_624.tsv'),sep='\t',row.names=FALSE,quote=TRUE)

## Exact original five-pathway composite models using frozen sample-level scores.
q<-d[complete.cases(d[,c('metabolic_intensity','cluster','proliferation')]),]
m0<-lm(metabolic_intensity~cluster,q);m1<-lm(metabolic_intensity~cluster+proliferation,q);c0<-summary(m0)$coef['cluster2',];c1<-summary(m1)$coef['cluster2',]
intres<-data.frame(model=c('Unadjusted','Proliferation-adjusted'),n=nrow(q),beta_C2_vs_C1=c(c0[1],c1[1]),SE=c(c0[2],c1[2]),p=c(c0[4],c1[4]),FDR=p.adjust(c(c0[4],c1[4]),'BH'),C1_mean=mean(q[q$cluster==1,'metabolic_intensity']),C2_mean=mean(q[q$cluster==2,'metabolic_intensity']),proliferation_beta=c(NA,summary(m1)$coef['proliferation','Estimate']))
write.table(intres,file.path(out,'five_pathway_metabolic_intensity_624.tsv'),sep='\t',row.names=FALSE,quote=TRUE)

## Old-versus-new comparison using frozen published numeric sources.
oldadj<-data.frame(outcome=c('lactate_axis','acidosis_composite'),old_n=645,old_beta=c(0.201079351714311,0.0820690052846843),old_FDR=c(1.39930811142847e-05,0.0180724322864827))
newadj<-adj[match(oldadj$outcome,adj$outcome),c('outcome','n','beta','FDR')];names(newadj)[2:4]<-c('new_n','new_beta','new_FDR');compadj<-merge(oldadj,newadj,by='outcome');write.table(compadj,file.path(out,'old_vs_new_acidification_adjusted.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
oldint<-data.frame(model=c('Unadjusted','Proliferation-adjusted'),old_n=647,old_beta=c(1.1544,0.8959),old_FDR=c(7.00e-65,3.54e-27));compint<-merge(oldint,intres[,c('model','n','beta_C2_vs_C1','FDR')],by='model');names(compint)[5:7]<-c('new_n','new_beta','new_FDR');write.table(compint,file.path(out,'old_vs_new_metabolic_intensity.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
writeLines(c('Patient-level sensitivity analysis using the frozen scores and exact original model formulas.','Acidification model: marker ~ acidosis_signed_score + proliferation_composite + hypoxia_score + microenv_proxy + CMS + project + stage_group.','Metabolic intensity model: metabolic_intensity ~ cluster + proliferation.','Scores were not recalculated or rescaled after deduplication; the fixed 624-patient table selected one frozen sample record per patient.'),file.path(out,'model_definition_and_scope.txt'))
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
