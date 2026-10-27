suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(png)
  library(grid)
  library(scales)
})

inp <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_results310_output/source_inputs")
out <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_results310_output")
dir.create(out,recursive=TRUE,showWarnings=FALSE)

ia <- fread(file.path(inp,"integration_audit.tsv"))
f8 <- fread(file.path(inp,"final8.tsv"))
texp <- fread(file.path(inp,"target_expression.tsv"))[target %in% c("WEE1","PLK1")]
bind <- fread(file.path(inp,"direct_binding_audit.tsv"))
struc <- fread(file.path(inp,"structure_audit.tsv"))
rmsd <- fread(file.path(inp,"redocking_rmsd.tsv"))
rdcon <- fread(file.path(inp,"redocking_conclusion.tsv"))
io <- fread(file.path(inp,"interaction_overlap.tsv"))

funnel <- data.table(
  step=factor(c("GDSC2 drugs","FDR < 0.05","Effect threshold","Priority mechanism","Evidence shortlist","Dual co-crystals"),
    levels=rev(c("GDSC2 drugs","FDR < 0.05","Effect threshold","Priority mechanism","Evidence shortlist","Dual co-crystals"))),
  n=c(198,164,80,39,8,1)
)

f8[,drug_label:=sub("_[0-9]+$","",drug)]
f8[,drug_label:=factor(drug_label,levels=drug_label[order(delta_median)])]
f8[,predicted_direction:=ifelse(delta_median<0,"C2 predicted more sensitive","C1 predicted more sensitive")]

score <- f8[,.(drug_label,refined_score,target_expression_score,PPI_score,structure_score,CRC_score,drug_effect_score)]
score[,candidate:=as.character(drug_label)]
score[,highlight:=candidate=="MK-1775"]
score[,candidate:=factor(candidate,levels=as.character(f8$drug_label)[order(f8$refined_score)] )]

selected_rmsd <- rmsd[(pdb_id=="5V5Y" & pose==1) | (pdb_id=="9D0P" & pose %in% c(1,2))]
selected_rmsd[,label:=factor(ifelse(pdb_id=="5V5Y","WEE1 pose 1",paste0("PLK1 pose ",pose)),
  levels=c("WEE1 pose 1","PLK1 pose 2","PLK1 pose 1"))]

contacts <- io[interaction_type=="ALL_TYPED_CONTACTS"]
contacts[,label:=factor(c("WEE1 pose 1","PLK1 pose 1","PLK1 pose 2"),levels=c("WEE1 pose 1","PLK1 pose 2","PLK1 pose 1"))]

# Archive panel-ready data and evidence audits.
fwrite(funnel,file.path(out,"Figure7A_drug_screening_funnel.tsv"),sep="\t")
fwrite(f8,file.path(out,"Figure7B_final8_predicted_sensitivity.tsv"),sep="\t")
fwrite(score,file.path(out,"Figure7C_final8_evidence_scores.tsv"),sep="\t")
fwrite(texp,file.path(out,"Figure7D_WEE1_PLK1_expression.tsv"),sep="\t")
fwrite(selected_rmsd,file.path(out,"Figure7E_selected_redocking_RMSD.tsv"),sep="\t")
fwrite(contacts,file.path(out,"Figure7F_interaction_retention.tsv"),sep="\t")
fwrite(bind,file.path(out,"direct_binding_and_cocrystal_audit.tsv"),sep="\t")
fwrite(struc,file.path(out,"WEE1_PLK1_structure_audit.tsv"),sep="\t")
fwrite(rdcon,file.path(out,"redocking_environment_and_conclusion.tsv"),sep="\t")

theme_fig <- theme_minimal(base_size=9.1)+theme(panel.grid.minor=element_blank(),
  plot.title=element_text(face="bold",size=10.2),plot.subtitle=element_text(size=7.6),
  axis.title=element_text(size=8.3),legend.position="bottom",legend.title=element_blank(),
  plot.margin=margin(6,7,6,7))

pA <- ggplot(funnel,aes(n,step,fill=n))+geom_col(width=.68)+geom_text(aes(label=n),hjust=-.3,size=2.7,fontface="bold")+
  scale_fill_gradient(low="#B9D3C2",high="#315B4C",guide="none")+scale_x_continuous(limits=c(0,218))+
  labs(title="A  Candidate-screening funnel",subtitle="Prediction and evidence prioritization; no experimental efficacy claim",x="Candidates retained",y=NULL)+theme_fig

pB <- ggplot(f8,aes(delta_median,drug_label,fill=predicted_direction))+geom_vline(xintercept=0,color="#777777")+
  geom_col(width=.68)+scale_fill_manual(values=c(`C2 predicted more sensitive`="#3E7CB1",`C1 predicted more sensitive`="#D77A61"))+
  labs(title="B  Eight shortlisted drug predictions",subtitle="C2-C1 predicted response; negative = predicted C2 sensitivity",x="Median predicted response difference",y=NULL)+theme_fig+
  theme(axis.text.y=element_text(size=7.2),legend.text=element_text(size=7.1))

pC <- ggplot(score,aes(refined_score,candidate,fill=highlight))+geom_col(width=.68)+geom_text(aes(label=refined_score),hjust=-.35,size=2.6)+
  scale_fill_manual(values=c(`TRUE`="#B54A4A",`FALSE`="#8AA6A3"),guide="none")+scale_x_continuous(limits=c(0,10.2))+
  labs(title="C  Integrated evidence priority",subtitle="Expression, PPI, structure, CRC evidence, and drug-effect scoring",x="Refined evidence score",y=NULL)+theme_fig+
  theme(axis.text.y=element_text(size=7.2))

pD <- ggplot(texp,aes(TCGA_delta_C2_C1,target,fill=target))+geom_vline(xintercept=0,color="#777777")+geom_col(width=.62)+
  geom_text(aes(label=paste0("FDR=",format(TCGA_FDR,digits=2,scientific=TRUE))),hjust=-.12,size=2.55)+
  scale_fill_manual(values=c(WEE1="#4C78A8",PLK1="#E07A5F"),guide="none")+scale_x_continuous(limits=c(0,29))+
  labs(title="D  Target expression in C2 versus C1",subtitle="TCGA TPM difference; both targets detected in all samples",x="C2-C1 median TPM",y=NULL)+theme_fig

pE <- ggplot(selected_rmsd,aes(heavy_atom_RMSD_A,label,fill=passes_RMSD_le_2A))+geom_vline(xintercept=2,linetype=2,color="#B54A4A")+
  geom_col(width=.64)+geom_text(aes(label=sprintf("%.3f A",heavy_atom_RMSD_A)),hjust=ifelse(selected_rmsd$heavy_atom_RMSD_A<2,-.15,1.1),color=ifelse(selected_rmsd$heavy_atom_RMSD_A<2,"#333333","white"),size=2.55)+
  scale_fill_manual(values=c(`TRUE`="#4C956C",`FALSE`="#B85C5C"),labels=c(`TRUE`="RMSD <= 2 A",`FALSE`="RMSD > 2 A"))+
  scale_x_continuous(limits=c(0,9.4))+labs(title="E  Co-crystal redocking",subtitle="Dashed line: prespecified 2-A heavy-atom RMSD threshold",x="Heavy-atom RMSD (A)",y=NULL)+theme_fig+
  theme(legend.text=element_text(size=7.1))

pF <- ggplot(contacts,aes(reference_retention_fraction,label,fill=target))+geom_col(width=.64)+
  geom_text(aes(label=percent(reference_retention_fraction,accuracy=.1)),hjust=-.2,size=2.6)+
  scale_fill_manual(values=c(WEE1="#4C78A8",PLK1="#E07A5F"),guide="none")+scale_x_continuous(limits=c(0,1.12),labels=percent_format())+
  labs(title="F  Crystal-contact retention",subtitle="All typed PLIP contacts retained in selected redocked poses",x="Reference-contact retention",y=NULL)+theme_fig

imgW <- readPNG(file.path(inp,"WEE1_overlay.png")); imgP <- readPNG(file.path(inp,"PLK1_overlay.png"))
pG <- ggplot()+annotation_custom(rasterGrob(imgW,interpolate=TRUE),-Inf,Inf,-Inf,Inf)+coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  labs(title="G  WEE1: crystal versus redocked pose")+theme_void()+theme(plot.title=element_text(face="bold",size=10.2),plot.margin=margin(6,7,6,7))
pH <- ggplot()+annotation_custom(rasterGrob(imgP,interpolate=TRUE),-Inf,Inf,-Inf,Inf)+coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  labs(title="H  PLK1: ranking limitation and native-like pose")+theme_void()+theme(plot.title=element_text(face="bold",size=10.2),plot.margin=margin(6,7,6,7))

fig7 <- (pA|pB)/(pC|pD)/(pE|pF)/(pG|pH)+plot_layout(heights=c(.9,.9,.9,1.6))+
  plot_annotation(caption="Exploratory computational evidence only. Predictions and redocking do not establish clinical efficacy, dual-target pharmacology, or causal vulnerability.\nNo unfinished MD result is included.",
    theme=theme(plot.caption=element_text(size=7.6,color="#8B2E2E",hjust=.5,face="italic",lineheight=.95),plot.margin=margin(4,4,10,4)))

ggsave(file.path(out,"Figure7_analysis_version.pdf"),fig7,width=210,height=300,units="mm",device=cairo_pdf)
ggsave(file.path(out,"Figure7_analysis_version_300dpi.png"),fig7,width=210,height=300,units="mm",dpi=300,bg="white")

audit <- data.table(section="3.10",metric=c("total_drugs","FDR_significant","effect_filtered","mechanism_priority","evidence_shortlist","dual_cocrystal_candidate",
  "MK1775_delta_C2_C1","MK1775_FDR","WEE1_redock_RMSD","PLK1_native_like_RMSD","WEE1_contact_retention","PLK1_native_like_contact_retention"),
  value=c(198,164,80,39,8,1,f8[grepl("MK-1775",drug),delta_median],f8[grepl("MK-1775",drug),FDR],
    selected_rmsd[pdb_id=="5V5Y",heavy_atom_RMSD_A],selected_rmsd[pdb_id=="9D0P" & pose==2,heavy_atom_RMSD_A],
    contacts[target=="WEE1",reference_retention_fraction],contacts[target=="PLK1" & grepl("pose2",comparison),reference_retention_fraction]))
fwrite(audit,file.path(out,"Results_3.10_frozen_number_audit.tsv"),sep="\t")

writeLines(c(
  "Figure 7. Exploratory drug-response prediction and structural compatibility analyses prioritize a testable MK-1775-WEE1/PLK1 hypothesis.",
  "(A) Sequential drug-screening and evidence-prioritization funnel.",
  "(B) Median C2-minus-C1 predicted response differences for the eight shortlisted drugs. These are oncoPredict estimates, not measured drug responses.",
  "(C) Integrated evidence scores combining target expression, PPI, structure availability, CRC relevance, and predicted drug effect.",
  "(D) TCGA WEE1 and PLK1 expression differences between C2 and C1.",
  "(E) Heavy-atom RMSD for selected WEE1 and PLK1 redocked poses. The PLK1 top-ranked pose failed the threshold, whereas pose 2 was native like.",
  "(F) Fraction of crystallographic PLIP contacts retained by the selected redocked poses.",
  "(G-H) Ligand-coordinate overlays for WEE1 and PLK1. Structural panels show pose compatibility only and do not demonstrate biochemical dual-target inhibition."
),file.path(out,"Figure7_legend_draft.md"))

capture.output(sessionInfo(),file=file.path(out,"R_sessionInfo.txt"))
cat("Figure 7 completed\n")
