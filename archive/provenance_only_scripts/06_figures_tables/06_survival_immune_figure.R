suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(survival)
  library(scales)
})

inp <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_results39/inputs")
out <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_results39_output")
dir.create(out,recursive=TRUE,showWarnings=FALSE)

# The historical survival summary contains literal `t separators; normalize before parsing.
surv_lines <- readLines(file.path(inp,"survival_summary_raw.tsv"),warn=FALSE)
surv_txt <- paste(gsub("`t","\t",surv_lines,fixed=TRUE),collapse="\n")
surv_all <- fread(text=surv_txt)
surv2 <- surv_all[k==2]
surv2[, endpoint := factor(endpoint,levels=c("OS","PFI","DFI","DSS"))]

clin <- fread(file.path(inp,"k2_survival_samples.tsv"),check.names=FALSE)
clin[, state := factor(paste0("C",cluster),levels=c("C1","C2"))]
os <- clin[!is.na(OS) & !is.na(OS.time) & OS.time>=0,.(sample_id,state,time_months=as.numeric(OS.time)/30.4375,status=as.numeric(OS))]
fit <- survfit(Surv(time_months,status)~state,data=os)
sfit <- summary(fit)
km <- data.table(time_months=sfit$time,survival=sfit$surv,state=sub("^state=","",sfit$strata))
km <- rbind(data.table(time_months=0,survival=1,state=c("C1","C2")),km)
km[,state:=factor(state,levels=c("C1","C2"))]

mcp <- fread(file.path(inp,"mcp.tsv")); xc <- fread(file.path(inp,"xcell.tsv"))
cross <- fread(file.path(inp,"robust.tsv"))
crosslong <- rbind(
  cross[,.(feature,method="MCP-counter",matched_label=feature,effect_C2_minus_C1=-MCP_delta,FDR=MCP_FDR)],
  cross[,.(feature,method="xCell",matched_label=method2,effect_C2_minus_C1=-xCell_delta,FDR=xCell_FDR)]
)
crosslong[,direction:=fifelse(effect_C2_minus_C1>0,"C2 higher","C1 higher")]
crosslong[,sig:=FDR<.05]

c7 <- fread(file.path(inp,"c7_theme.tsv"))
c7[,effect_C2_minus_C1 := -c7_median_delta]
c7[,theme_label:=factor(theme,levels=rev(c("CD8_T","T_cell","B_cell","NK","Myeloid","Interferon","Immune_suppression")),
  labels=rev(c("CD8 T","T cell","B cell","NK","Myeloid","Interferon","Immune suppression")))]

chk <- fread(file.path(inp,"checkpoint.tsv"))
chk[,`:=`(effect_C2_minus_C1=-delta_C1_minus_C2,signed_significance=sign(-delta_C1_minus_C2)*(-log10(pmax(FDR,1e-300))))]
chk[,gene:=factor(gene,levels=gene[order(signed_significance)])]

# Panel-ready archive.
fwrite(surv2,file.path(out,"Figure6A_K2_four_endpoint_survival.tsv"),sep="\t")
fwrite(os,file.path(out,"Figure6B_OS_sample_level.tsv.gz"),sep="\t")
fwrite(km,file.path(out,"Figure6B_OS_KM_coordinates.tsv"),sep="\t")
fwrite(crosslong,file.path(out,"Figure6C_MCP_xCell_direction_matrix.tsv"),sep="\t")
fwrite(c7,file.path(out,"Figure6D_C7_theme_summary.tsv"),sep="\t")
fwrite(chk,file.path(out,"Figure6E_checkpoint_effector_differential.tsv"),sep="\t")

theme_fig <- theme_minimal(base_size=9.3)+theme(panel.grid.minor=element_blank(),
  plot.title=element_text(face="bold",size=10.4),plot.subtitle=element_text(size=7.8),
  axis.title=element_text(size=8.5),legend.position="bottom",legend.title=element_blank(),
  plot.margin=margin(6,7,6,7))

# A: all four endpoints are non-significant.
surv_x_labels <- setNames(paste0(as.character(surv2$endpoint),"\nn=",surv2$n,", e=",surv2$events),as.character(surv2$endpoint))
pA <- ggplot(surv2,aes(endpoint,-log10(logrank_p),fill=endpoint))+geom_hline(yintercept=-log10(.05),linetype=2,color="#B54A4A")+
  geom_col(width=.65)+geom_text(aes(label=paste0("P=",sprintf("%.3f",logrank_p))),vjust=-.4,size=2.5)+
  scale_fill_manual(values=c("#4C78A8","#72A0C1","#A6A6A6","#7A7A7A"),guide="none")+
  scale_x_discrete(labels=surv_x_labels)+scale_y_continuous(limits=c(0,1.65),breaks=c(0,.5,1,-log10(.05),1.5))+
  labs(title="A  Four survival endpoints",subtitle="K=2 log-rank tests; dashed line denotes P=0.05",x=NULL,y="-log10(P)")+theme_fig

# B: representative OS curve.
pB <- ggplot(km,aes(time_months,survival,color=state))+geom_step(linewidth=.75)+
  scale_color_manual(values=c(C1="#2F6690",C2="#D77A61"))+scale_y_continuous(labels=percent_format(accuracy=1),limits=c(0,1.02))+
  coord_cartesian(xlim=c(0,150))+labs(title="B  Overall survival",subtitle="Representative Kaplan-Meier curve; log-rank P=0.198 (n=643)",x="Months",y="Survival probability")+theme_fig

# C: matched MCP-counter/xCell evidence.
feature_order <- rev(c("CD8 T cells","B lineage","T cells","NK cells","Neutrophils","Endothelial cells","Fibroblasts"))
crosslong[,feature:=factor(feature,levels=feature_order)]
pC <- ggplot(crosslong,aes(method,feature,fill=direction))+geom_tile(color="white",linewidth=.7)+
  geom_text(aes(label=ifelse(sig,"*","")),size=4,fontface="bold")+
  scale_fill_manual(values=c(`C2 higher`="#D77A61",`C1 higher`="#3E7CB1"))+
  labs(title="C  MCP-counter/xCell cross-validation",subtitle="Direction is C2-C1; * FDR<0.05. No mapped feature met the robust rule",x=NULL,y=NULL)+theme_fig+
  theme(panel.grid=element_blank(),axis.text.x=element_text(angle=15,hjust=1))

# D: C7 theme-level evidence.
pD <- ggplot(c7,aes(effect_C2_minus_C1,theme_label))+geom_vline(xintercept=0,linetype=2,color="#777777")+
  geom_col(aes(fill=c7_fdr05>0),width=.65)+geom_text(aes(label=paste0(c7_fdr05,"/",c7_n," FDR<0.05")),hjust=-.15,size=2.5)+
  scale_fill_manual(values=c(`TRUE`="#6A9A8B",`FALSE`="#B9B9B9"),guide="none")+
  scale_x_continuous(limits=c(0,.34))+labs(title="D  C7 immune themes",subtitle="Median C2-C1 GSVA effect; not a cell-fraction estimate",x="Median GSVA difference",y=NULL)+theme_fig

# E: checkpoint and effector genes shown on a signed significance scale.
pE <- ggplot(chk,aes(signed_significance,gene,fill=FDR<.05))+geom_vline(xintercept=0,color="#777777")+geom_col(width=.68)+
  scale_fill_manual(values=c(`TRUE`="#C95D63",`FALSE`="#B9B9B9"),labels=c(`TRUE`="FDR < 0.05",`FALSE`="Not significant"))+
  labs(title="E  Immune checkpoint and effector genes",subtitle="Signed -log10(FDR); positive values indicate higher expression in C2",x="Signed -log10(FDR)",y=NULL)+theme_fig+
  theme(axis.text.y=element_text(size=7.1),legend.text=element_text(size=7.2))

# F: explicit clinical boundary.
pF <- ggplot()+annotate("rect",xmin=.04,xmax=.96,ymin=.08,ymax=.92,fill="#F3F5F6",color="#64717A",linewidth=.65)+
  annotate("text",x=.5,y=.72,label="Biological-state evidence",fontface="bold",size=3.6)+
  annotate("text",x=.5,y=.53,label="Method-dependent immune differences\nwithout a three-method robust cell type",size=3.05,lineheight=.95)+
  annotate("text",x=.5,y=.28,label="No stable prognostic separation\nacross OS, PFI, DFI, or DSS",size=3.05,lineheight=.95,color="#A23B3B")+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))+labs(title="F  Clinical interpretation boundary")+theme_void()+
  theme(plot.title=element_text(face="bold",size=10.4),plot.margin=margin(6,7,6,7))

fig6 <- (pA|pB)/(pC|pD)/(pE|pF)+plot_layout(heights=c(1,1,1.12))
ggsave(file.path(out,"Figure6_analysis_version.pdf"),fig6,width=205,height=225,units="mm",device=cairo_pdf)
ggsave(file.path(out,"Figure6_analysis_version_300dpi.png"),fig6,width=205,height=225,units="mm",dpi=300,bg="white")

audit <- data.table(section=c(rep("survival",4),"immune","immune","immune","immune"),
  metric=c(as.character(surv2$endpoint),"MCP_FDR05_n","xCell_FDR05_n","mapped_robust_n","checkpoint_effector_FDR05_n"),
  value=c(surv2$logrank_p,sum(mcp$FDR<.05),sum(xc$FDR<.05),sum(cross$robust),sum(chk$FDR<.05)))
fwrite(audit,file.path(out,"Results_3.9_frozen_number_audit.tsv"),sep="\t")

writeLines(c(
  "Figure 6. Survival and immune analyses define the clinical boundary of the frozen expression state.",
  "(A) Log-rank results for OS, PFI, DFI, and DSS using the K=2 state labels. The horizontal line marks nominal P=0.05.",
  "(B) Representative overall-survival Kaplan-Meier curves for C1 and C2.",
  "(C) Direction and FDR status for seven manually mapped MCP-counter/xCell feature pairs. Red indicates higher values in C2 and blue indicates higher values in C1.",
  "(D) Median C2-minus-C1 GSVA differences for prespecified C7 immune themes, with the number of significant signatures in each theme.",
  "(E) Immune-checkpoint and effector-gene differences on a signed -log10(FDR) scale. Positive values indicate higher expression in C2.",
  "(F) Interpretation boundary: immune findings were method dependent, and no stable survival separation was identified."
),file.path(out,"Figure6_legend_draft.md"))

capture.output(sessionInfo(),file=file.path(out,"R_sessionInfo.txt"))
cat("Figure 6 completed\n")
