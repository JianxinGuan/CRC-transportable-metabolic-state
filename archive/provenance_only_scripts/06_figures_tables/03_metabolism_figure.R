options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(scales)})
set.seed(20260830)
inp <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_results65_gse39582_bom_corrected/inputs"); out <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_results65_gse39582_bom_corrected")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

pw <- fread(file.path(inp,"pathways.tsv")); inten <- fread(file.path(inp,"intensity.tsv")); isamp <- fread(file.path(inp,"intensity_samples.tsv.gz"))
bal <- fread(file.path(inp,"balance.tsv")); cors <- fread(file.path(inp,"acid_spearman.tsv")); models <- fread(file.path(inp,"acid_models.tsv"))
k2 <- fread(file.path(inp,"acid_k2.tsv")); ext <- fread(file.path(inp,"acid_external.tsv")); asamp <- fread(file.path(inp,"acid_samples.tsv.gz")); ov <- fread(file.path(inp,"overlap.tsv"))

cohort_levels <- c("TCGA_COADREAD","GSE39582","GSE17536")
cohort_labels <- c(TCGA_COADREAD="TCGA-COAD/READ",GSE39582="GSE39582",GSE17536="GSE17536")
path_levels <- c("HALLMARK_GLYCOLYSIS","HALLMARK_MTORC1_SIGNALING","HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY","HALLMARK_FATTY_ACID_METABOLISM")
path_labels <- c(HALLMARK_GLYCOLYSIS="Glycolysis",HALLMARK_MTORC1_SIGNALING="mTORC1",HALLMARK_OXIDATIVE_PHOSPHORYLATION="OXPHOS",HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY="ROS pathway",HALLMARK_FATTY_ACID_METABOLISM="Fatty acid metabolism")

# Figure 3 data: all effects use C2 minus C1.
p5 <- pw[pathway %chin% path_levels]
p5[, `:=`(effect_C2_minus_C1=-delta_C1_minus_C2, cohort=factor(cohort,cohort_levels), pathway=factor(pathway,path_levels))]
p5[, sig := FDR_within_cohort < .05]
fwrite(p5,file.path(out,"Figure3_five_pathway_effects.tsv"),sep="\t")

intlong <- melt(inten[,.(cohort,unadjusted_beta_C2_vs_C1,unadjusted_FDR,adjusted_beta_C2_vs_C1,adjusted_FDR)],id.vars="cohort",
                measure.vars=list(beta=c("unadjusted_beta_C2_vs_C1","adjusted_beta_C2_vs_C1"),FDR=c("unadjusted_FDR","adjusted_FDR")),
                variable.name="model")
intlong[,model:=factor(model,levels=1:2,labels=c("Unadjusted","Proliferation-adjusted"))]
intlong[,cohort:=factor(cohort,cohort_levels)]
fwrite(intlong,file.path(out,"Figure3_metabolic_intensity_adjustment.tsv"),sep="\t")

isamp[,`:=`(cohort=factor(cohort,cohort_levels),cluster=factor(cluster,levels=c(1,2),labels=c("C1","C2")))]
fwrite(isamp,file.path(out,"Figure3_metabolic_intensity_sample_level.tsv.gz"),sep="\t")

bplot <- bal[outcome=="GLY_minus_OXPHOS_balance"]
bplot[,`:=`(cohort=factor(cohort,cohort_levels),model=factor("Unadjusted"),beta=unadjusted_beta_C2_vs_C1,FDR=NA_real_)]
bplot2 <- bal[outcome=="GLY_minus_OXPHOS_balance",.(cohort=factor(cohort,cohort_levels),model=factor("Proliferation-adjusted"),beta=adjusted_beta_C2_vs_C1,FDR=adjusted_FDR)]
bplot <- rbind(bplot[,.(cohort,model,beta,FDR)],bplot2); bplot[,model:=factor(as.character(model),levels=c("Unadjusted","Proliferation-adjusted"))]
fwrite(bplot,file.path(out,"Figure3_glycolysis_OXPHOS_balance.tsv"),sep="\t")

# Figure 4 data.
marker_labels <- c(carbonic_anhydrase="Carbonic anhydrases",acid_extrusion_transport="Acid/proton transport",lactate_axis="Lactate axis",ca9_ca12="CA9/CA12",acidosis_composite="Acidification composite",CA9="CA9",CA12="CA12")
cors[,label:=factor(marker_labels[outcome],levels=rev(marker_labels[c("carbonic_anhydrase","acid_extrusion_transport","lactate_axis","ca9_ca12","acidosis_composite","CA9","CA12")]))]
k2[,label:=factor(marker_labels[outcome],levels=rev(marker_labels[c("carbonic_anhydrase","acid_extrusion_transport","lactate_axis","ca9_ca12","acidosis_composite","CA9","CA12")]))]
fwrite(cors,file.path(out,"Figure4_TCGA_unadjusted_correlations.tsv"),sep="\t"); fwrite(k2,file.path(out,"Figure4_TCGA_K2_marker_differences.tsv"),sep="\t")

forest <- models[outcome %chin% c("acidosis_composite","lactate_axis","acid_extrusion_transport","carbonic_anhydrase","ca9_ca12")]
forest[,`:=`(label=factor(marker_labels[outcome],levels=rev(marker_labels[c("carbonic_anhydrase","acid_extrusion_transport","lactate_axis","ca9_ca12","acidosis_composite")])),analysis=factor(analysis,levels=c("all","exclude13"),labels=c("All evaluable","Exclude 13")))]
fwrite(forest,file.path(out,"Figure4_adjusted_models.tsv"),sep="\t")

extplot <- ext[outcome %chin% c("acidosis_composite","lactate_axis")]
extplot[,`:=`(label=factor(marker_labels[outcome],levels=c("Acidification composite","Lactate axis")),cohort=factor(cohort,levels=c("GSE39582","GSE17536")))]
fwrite(extplot,file.path(out,"Figure4_external_marker_effects.tsv"),sep="\t")

overlap_plot <- ov[,.(signature,n_predefined,n_detected,n_overlap68)]
overlap_plot[,label:=factor(c("Carbonic anhydrases","Acid/proton transport","Lactate axis","CA9/CA12"),levels=rev(c("Carbonic anhydrases","Acid/proton transport","Lactate axis","CA9/CA12")))]
fwrite(overlap_plot,file.path(out,"Figure4_marker_overlap_audit.tsv"),sep="\t")

# Frozen audit.
audit <- rbind(
  inten[,.(section="metabolic_intensity",cohort,metric=c("unadjusted_beta","adjusted_beta"),value=c(unadjusted_beta_C2_vs_C1,adjusted_beta_C2_vs_C1),FDR=c(unadjusted_FDR,adjusted_FDR)),by=cohort],
  models[outcome=="acidosis_composite",.(section="patient_acidification",cohort=analysis,metric="continuous_score_adjusted_beta",value=beta,FDR)],
  ext[outcome %chin% c("acidosis_composite","lactate_axis"),.(section="external_markers",cohort,metric=paste0(outcome,"_delta_C2_minus_C1"),value=delta_C2_minus_C1,FDR)]
,fill=TRUE)
fwrite(audit,file.path(out,"Results_3.5_3.6_frozen_number_audit.tsv"),sep="\t")

theme_fig <- theme_minimal(base_size=9.5)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=10.5),plot.subtitle=element_text(size=8),axis.title=element_text(size=8.5),legend.position="bottom",legend.title=element_blank(),plot.margin=margin(6,6,6,6))
cols <- c(C1="#2678A5",C2="#D65F5F")

# Figure 3.
p3A <- ggplot(p5,aes(cohort,pathway,fill=effect_C2_minus_C1))+geom_tile(color="white",linewidth=.7)+
  geom_text(aes(label=paste0(sprintf("%.2f",effect_C2_minus_C1),ifelse(sig,"*",""))),size=3)+
  scale_fill_gradient2(low="#3B6FB6",mid="white",high="#C94C4C",midpoint=0)+
  scale_x_discrete(labels=cohort_labels)+scale_y_discrete(labels=path_labels)+
  labs(title="A  Five-pathway direction map",subtitle="Effect = C2 - C1; * within-cohort FDR < 0.05",x=NULL,y=NULL,fill="Effect")+theme_fig+theme(axis.text.x=element_text(angle=20,hjust=1))
p3B <- ggplot(p5,aes(pathway,effect_C2_minus_C1,color=cohort,shape=sig))+geom_hline(yintercept=0,color="#777777",linetype=2)+geom_point(size=2.8)+
  scale_x_discrete(labels=path_labels)+scale_color_manual(values=c("#333333","#2A9D8F","#E07A5F"),labels=c("TCGA","GSE39582","GSE17536"))+scale_shape_manual(values=c(`TRUE`=16,`FALSE`=1))+
  labs(title="B  Cross-cohort pathway effects",subtitle="Filled points meet within-cohort FDR < 0.05",x=NULL,y="C2 - C1 mean GSVA score",color=NULL,shape=NULL)+guides(shape=guide_legend(order=1),color=guide_legend(order=2))+theme_fig+theme(axis.text.x=element_text(angle=25,hjust=1),legend.box="vertical",legend.text=element_text(size=7.4),legend.key.width=unit(3.5,"mm"))
p3C <- ggplot(isamp,aes(cluster,metabolic_intensity,fill=cluster))+geom_violin(trim=FALSE,alpha=.45,color=NA)+geom_boxplot(width=.18,outlier.shape=NA,linewidth=.35)+facet_wrap(~cohort,labeller=labeller(cohort=c(TCGA_COADREAD="TCGA",GSE39582="GSE39582",GSE17536="GSE17536")),scales="free_y")+
  scale_fill_manual(values=cols,guide="none")+labs(title="C  Five-pathway metabolic intensity",subtitle="Unadjusted sample-level distributions",x=NULL,y="Composite intensity")+theme_fig+theme(strip.text=element_text(size=7.3))
p3D <- ggplot(intlong,aes(cohort,beta,fill=model))+geom_hline(yintercept=0,linetype=2,color="#666666")+geom_col(position=position_dodge(.72),width=.65)+geom_text(aes(label=ifelse(FDR<.05,"*","ns")),position=position_dodge(.72),vjust=ifelse(intlong$beta>=0,-.5,1.3),size=2.7)+
  scale_x_discrete(labels=cohort_labels)+scale_fill_manual(values=c("#7896AD","#D07B5F"))+
  labs(title="D  Proliferation-adjustment sensitivity",subtitle="Regression beta for C2 vs C1; * FDR < 0.05",x=NULL,y="Beta (C2 vs C1)")+theme_fig+theme(axis.text.x=element_text(angle=18,hjust=1))
p3E <- ggplot(bplot,aes(cohort,beta,fill=model))+geom_hline(yintercept=0,linetype=2,color="#666666")+geom_col(position=position_dodge(.72),width=.65)+
  scale_x_discrete(labels=cohort_labels)+scale_fill_manual(values=c("#7896AD","#D07B5F"))+
  labs(title="E  Glycolysis-OXPHOS balance",subtitle="No reproducible cross-cohort switching pattern",x=NULL,y="Beta (C2 vs C1)")+theme_fig+theme(axis.text.x=element_text(angle=18,hjust=1))
fig3 <- (p3A|p3B)/(p3C|p3D)/p3E+plot_layout(heights=c(1,1,1))
ggsave(file.path(out,"Figure3_analysis_version.pdf"),fig3,width=200,height=230,units="mm",device=cairo_pdf)
ggsave(file.path(out,"Figure3_analysis_version_300dpi.png"),fig3,width=200,height=230,units="mm",dpi=300,bg="white")

# Figure 4.
p4A <- ggplot(overlap_plot,aes(n_detected,label))+geom_col(fill="#5B8E7D",width=.65)+geom_text(aes(label=paste0(n_detected," detected\n0 overlap")),hjust=-.08,size=2.8)+scale_x_continuous(expand=expansion(mult=c(0,.28)))+
  labs(title="A  Independent marker-set audit",subtitle="All markers detected; zero overlap with the frozen 68 genes",x="Genes",y=NULL)+theme_fig
p4B <- ggplot(k2,aes(delta_C2_minus_C1,label,fill=FDR<.05))+geom_vline(xintercept=0,linetype=2,color="#777777")+geom_col(width=.65)+scale_fill_manual(values=c(`TRUE`="#C95D63",`FALSE`="#B9B9B9"),guide="none")+
  labs(title="B  TCGA categorical differences",subtitle="Unadjusted effect = C2 - C1",x="Marker-score difference",y=NULL)+theme_fig
p4C <- ggplot(cors,aes(rho,label,fill=FDR<.05))+geom_vline(xintercept=0,linetype=2,color="#777777")+geom_col(width=.65)+scale_fill_manual(values=c(`TRUE`="#557A95",`FALSE`="#B9B9B9"),guide="none")+
  labs(title="C  TCGA unadjusted correlations",subtitle="Spearman correlation with continuous state score",x="Spearman rho",y=NULL)+theme_fig
p4D <- ggplot(forest,aes(beta,label,color=analysis))+geom_vline(xintercept=0,linetype=2,color="#777777")+geom_errorbarh(aes(xmin=CI_low,xmax=CI_high),height=.18,position=position_dodge(.45))+geom_point(size=2.3,position=position_dodge(.45))+
  scale_color_manual(values=c("#2F6690","#D77A61"))+labs(title="D  Fully adjusted TCGA models",subtitle="Beta per 1-unit state-score increase",x="Adjusted beta (95% CI)",y=NULL,color=NULL)+theme_fig
p4E <- ggplot(extplot,aes(delta_C2_minus_C1,label,fill=cohort))+geom_vline(xintercept=0,linetype=2,color="#777777")+geom_col(position=position_dodge(.72),width=.65)+geom_text(aes(label=ifelse(FDR<.05,"*","ns")),position=position_dodge(.72),hjust=-.4,size=2.8)+
  scale_fill_manual(values=c("#2A9D8F","#E07A5F"))+labs(title="E  External marker differences",subtitle="Effect = C2 - C1; * FDR < 0.05",x="Marker-score difference",y=NULL)+theme_fig
p4F <- ggplot()+annotate("rect",xmin=.04,xmax=.96,ymin=.10,ymax=.90,fill="#F3F5F6",color="#64717A",linewidth=.6)+
  annotate("text",x=.5,y=.69,label="Patient-side transcriptomic\nassociation",fontface="bold",size=3.4,lineheight=.95)+
  annotate("text",x=.5,y=.45,label="Supports an acidosis-adaptation-\nrelated label",size=2.9,lineheight=.95)+
  annotate("text",x=.5,y=.23,label="Does not measure tumor pH\nor establish causality",size=2.8,lineheight=.95,color="#A23B3B")+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))+labs(title="F  Interpretation boundary")+theme_void()+theme(plot.title=element_text(face="bold",size=10.5),plot.margin=margin(6,6,6,6))
fig4 <- (p4A|p4B)/(p4C|p4D)/(p4E|p4F)+plot_layout(heights=c(1,1.05,.9))
ggsave(file.path(out,"Figure4_analysis_version.pdf"),fig4,width=200,height=225,units="mm",device=cairo_pdf)
ggsave(file.path(out,"Figure4_analysis_version_300dpi.png"),fig4,width=200,height=225,units="mm",dpi=300,bg="white")

writeLines(c("Figure 3. The frozen state is associated with coordinated metabolic activity rather than a reproducible glycolysis-to-OXPHOS switch.",
"(A) C2-minus-C1 differences for five prespecified Hallmark pathways across TCGA-COAD/READ, GSE39582 and GSE17536. Asterisks indicate within-cohort BH FDR<0.05.",
"(B) Cross-cohort comparison of the same pathway effects. Positive values indicate higher activity in C2.",
"(C) Sample-level five-pathway composite intensity by projected state end.",
"(D) Unadjusted and proliferation-adjusted C2-versus-C1 coefficients. Adjustment retained a positive effect in TCGA, attenuated the GSE39582 effect to null and reversed the GSE17536 effect.",
"(E) Glycolysis-minus-OXPHOS balance coefficients, which did not show a stable cross-cohort switching pattern."),file.path(out,"Figure3_legend_draft.md"))
writeLines(c("Figure 4. Independent acid-handling and lactate-related transcriptomic markers provide qualified patient-side support for the state label.",
"(A) Predefined marker sets were fully detected and had zero overlap with the frozen 68-gene set.",
"(B) Unadjusted TCGA C2-minus-C1 marker-score differences.",
"(C) Unadjusted Spearman correlations between marker scores and the continuous state score.",
"(D) Fully adjusted TCGA linear models in all evaluable samples and after excluding 13 source/vial-associated samples. Models adjusted for proliferation, Hallmark hypoxia, xCell microenvironment proxy, CMS, COAD/READ and stage.",
"(E) External C2-minus-C1 differences for the acidification composite and lactate axis.",
"(F) Interpretation boundary. These analyses support transcriptomic association but do not directly measure intratumoral pH or establish causal acidosis adaptation."),file.path(out,"Figure4_legend_draft.md"))
capture.output(sessionInfo(),file=file.path(out,"R_sessionInfo.txt")); cat("Figures 3 and 4 completed\n")
