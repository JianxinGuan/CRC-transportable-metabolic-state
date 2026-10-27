suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

inp <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_results37_38/inputs")
out <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_results37_38_output")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

con <- fread(file.path(inp,"concordance.tsv"))
boot <- fread(file.path(inp,"bootstrap_2000.tsv.gz"))
bootsum <- fread(file.path(inp,"bootstrap_summary.tsv"))
bootind <- fread(file.path(inp,"individual_bootstrap.tsv"))
sil <- fread(file.path(inp,"silhouette.tsv"))[analysis != "Residual68_K2_exclusion13"]
resk <- fread(file.path(inp,"residual_k.tsv"))
cmsn <- fread(file.path(inp,"cms_counts.tsv"))
cmsg <- fread(file.path(inp,"cms_global.tsv"))
cmss <- fread(file.path(inp,"cms_sensitivity.tsv"))
within <- fread(file.path(inp,"within_cms.tsv"))[outcome == "metabolic_intensity"]
inc <- fread(file.path(inp,"incremental.tsv"))[outcome == "metabolic_intensity"]
inter <- fread(file.path(inp,"interaction.tsv"))[outcome == "metabolic_intensity"]

analysis_labels <- c(
  `47_gene_deletion_sensitivity`="47-gene deletion",
  E2F_G2M_full="E2F+G2M full",
  E2F_G2M_excluding_21_overlap="E2F+G2M exclusive",
  frozen68_residualized_on_E2F_G2M="68-gene residualized"
)
con[, label := factor(analysis_labels[analysis], levels=unname(analysis_labels))]

bootlong <- melt(boot, id.vars="replicate",
  measure.vars=c("delta_47_minus_full","delta_47_minus_exclusive"),
  variable.name="contrast", value.name="delta_ARI")
bootlong[, label := factor(contrast,
  levels=c("delta_47_minus_full","delta_47_minus_exclusive"),
  labels=c("47-gene minus full E2F+G2M","47-gene minus exclusive E2F+G2M"))]

sil_labels <- c(Frozen68_K2="Frozen 68 genes",Deletion47_K2="47 genes after deletion",
  E2F_G2M_full_K2="Full E2F+G2M",E2F_G2M_exclusive_K2="Exclusive E2F+G2M")
sil[, label := factor(sil_labels[analysis], levels=unname(sil_labels))]

cmslong <- melt(cmsn,id.vars="cluster",variable.name="CMS",value.name="n")
cmslong[, state := factor(paste0("C",cluster),levels=c("C1","C2"))]
cmslong[, proportion := n/sum(n),by=state]

within[, CMS := factor(CMS,levels=c("CMS1","CMS2","CMS3","CMS4"))]
inc[, analysis_label := factor(analysis,levels=c("all","exclude13"),labels=c("All evaluable","Exclude 13"))]

# Archive panel-ready tables.
fwrite(con,file.path(out,"Figure5A_proliferation_concordance.tsv"),sep="\t")
fwrite(bootlong,file.path(out,"Figure5B_paired_bootstrap_ARI_differences.tsv.gz"),sep="\t")
fwrite(bootsum,file.path(out,"Figure5B_bootstrap_summary.tsv"),sep="\t")
fwrite(sil,file.path(out,"Figure5C_silhouette_comparison.tsv"),sep="\t")
fwrite(resk,file.path(out,"Figure5D_residualized_K2_K6_stability.tsv"),sep="\t")
fwrite(cmslong,file.path(out,"Figure5E_CMS_contingency_long.tsv"),sep="\t")
fwrite(cmsg,file.path(out,"Figure5E_CMS_global_test.tsv"),sep="\t")
fwrite(cmss,file.path(out,"Figure5E_exclusion13_CMS_sensitivity.tsv"),sep="\t")
fwrite(within,file.path(out,"Figure5F_within_CMS_metabolic_effects.tsv"),sep="\t")
fwrite(inc,file.path(out,"Figure5G_CMS_incremental_models.tsv"),sep="\t")
fwrite(inter,file.path(out,"Figure5G_state_by_CMS_interaction.tsv"),sep="\t")

theme_fig <- theme_minimal(base_size=9.2)+theme(
  panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=10.3),
  plot.subtitle=element_text(size=7.8),axis.title=element_text(size=8.4),
  legend.position="bottom",legend.title=element_blank(),plot.margin=margin(6,7,6,7))

# A: agreement with the frozen K=2 labels.
conlong <- melt(con,id.vars=c("analysis","label"),measure.vars=c("ARI","NMI","agreement"),
  variable.name="metric",value.name="value")
conlong[, metric := factor(metric,levels=c("ARI","NMI","agreement"),labels=c("ARI","NMI","Agreement"))]
pA <- ggplot(conlong,aes(label,value,fill=metric))+geom_col(position=position_dodge(.72),width=.66)+
  geom_text(aes(label=sprintf("%.2f",value)),position=position_dodge(.72),vjust=-.35,size=2.35)+
  scale_fill_manual(values=c("#2F6690","#79A7A0","#D77A61"))+scale_y_continuous(limits=c(0,1.08),breaks=seq(0,1,.25))+
  labs(title="A  Proliferation-control concordance",subtitle="Agreement with frozen 68-gene K=2 labels (n=647)",x=NULL,y="Concordance")+
  theme_fig+theme(axis.text.x=element_text(angle=18,hjust=1),legend.key.width=unit(4,"mm"))

# B: paired bootstrap difference.
pB <- ggplot(bootlong,aes(delta_ARI,fill=label,color=label))+geom_density(alpha=.28,linewidth=.7)+
  geom_vline(xintercept=0,linetype=2,color="#666666")+
  geom_vline(data=bootsum,aes(xintercept=delta_median),color=c("#2F6690","#D77A61"),linewidth=.7)+
  scale_fill_manual(values=c("#2F6690","#D77A61"))+scale_color_manual(values=c("#2F6690","#D77A61"))+
  labs(title="B  Paired bootstrap ARI advantage",subtitle="2,000 paired sample bootstraps; both 95% CIs exclude zero",x="ARI difference",y="Density")+
  theme_fig+theme(legend.text=element_text(size=7.1))

# C: silhouette controls.
pC <- ggplot(sil,aes(label,silhouette_mean,fill=label))+geom_col(width=.65)+
  geom_text(aes(label=sprintf("%.3f",silhouette_mean)),vjust=-.45,size=2.7)+
  scale_fill_manual(values=c("#3E7CB1","#78A1BB","#D08C60","#B56B45"),guide="none")+
  scale_y_continuous(limits=c(0,.105),breaks=c(0,.025,.05,.075,.1))+
  labs(title="C  Geometric separation remains weak",subtitle="All mean silhouette widths are below 0.10",x=NULL,y="Mean silhouette")+
  theme_fig+theme(axis.text.x=element_text(angle=18,hjust=1))

# D: residualized K report.
pD <- ggplot(resk,aes(K,PAC_0.1_0.9))+geom_line(color="#6C757D",linewidth=.7)+geom_point(aes(size=min_cluster_n),color="#B54A4A")+
  geom_text(aes(label=paste0("min n=",min_cluster_n),hjust=ifelse(K==2,0,ifelse(K==6,1,.5))),vjust=-.9,size=2.5)+
  scale_x_continuous(breaks=2:6,expand=expansion(mult=c(.06,.06)))+scale_y_continuous(limits=c(0,.72))+
  scale_size_continuous(range=c(2.3,2.3),guide="none")+
  labs(title="D  Residualized K=2-6 stability",subtitle="A 13-sample minimum cluster recurred at every K",x="K",y="PAC (0.1-0.9)")+theme_fig

# E: external CMS labels, n=503.
pE <- ggplot(cmslong,aes(state,proportion,fill=CMS))+geom_col(width=.68,color="white",linewidth=.25)+
  geom_text(aes(label=ifelse(proportion>.08,paste0(CMS,"\n",n),ifelse(proportion>.05,as.character(n),""))),position=position_stack(vjust=.5),size=2.35,color="white",fontface="bold")+
  scale_fill_manual(values=c(CMS1="#6A4C93",CMS2="#1982C4",CMS3="#FFCA3A",CMS4="#8AC926"))+
  scale_y_continuous(labels=percent_format(accuracy=1),expand=expansion(mult=c(0,.04)))+
  labs(title="E  Association with external CMS labels",subtitle="n=503; Fisher P=1.77 x 10^-5; Cramer's V=0.219",x=NULL,y="CMS composition")+theme_fig

# F: within predicted CMS metabolic effects, n=645.
pF <- ggplot(within,aes(hedges_g,CMS))+geom_vline(xintercept=0,linetype=2,color="#777777")+
  geom_errorbar(aes(xmin=hedges_g_boot_CI_low,xmax=hedges_g_boot_CI_high),orientation="y",width=.18,color="#35607D")+
  geom_point(size=2.5,color="#C65D3B")+geom_text(aes(label=paste0("n=",n)),hjust=-.45,size=2.5)+
  scale_x_continuous(limits=c(0,3.9))+labs(title="F  Metabolic gradient within each CMS",subtitle="Predicted CMS n=645; Hedges' g (C2-C1), bootstrap 95% CI",x="Hedges' g",y=NULL)+theme_fig

# G: incremental information beyond CMS and covariates.
inclong <- melt(inc,id.vars=c("analysis","analysis_label","n","partial_R2","LRT_p"),
  measure.vars=c("adjusted_R2_base","adjusted_R2_full"),variable.name="model",value.name="adjusted_R2")
inclong[, model := factor(model,levels=c("adjusted_R2_base","adjusted_R2_full"),labels=c("CMS + covariates","+ continuous state score"))]
pG <- ggplot(inclong,aes(analysis_label,adjusted_R2,fill=model))+geom_col(position=position_dodge(.72),width=.65)+
  geom_text(aes(label=sprintf("%.3f",adjusted_R2)),position=position_dodge(.72),vjust=-.45,size=2.5)+
  scale_fill_manual(values=c("#9AA6B2","#3B7EA1"))+scale_y_continuous(limits=c(0,.64))+
  labs(title="G  Incremental metabolic information",subtitle="Adjusted R2 gain: 0.126 (n=645) and 0.101 after excluding 13 (n=632)",x=NULL,y="Adjusted R2")+
  theme_fig+theme(axis.text.x=element_text(angle=12,hjust=1),legend.text=element_text(size=7.2))

fig5 <- (pA|pB)/(pC|pD)/(pE|pF)/pG + plot_layout(heights=c(1,1,1.03,.9))
ggsave(file.path(out,"Figure5_analysis_version.pdf"),fig5,width=205,height=270,units="mm",device=cairo_pdf)
ggsave(file.path(out,"Figure5_analysis_version_300dpi.png"),fig5,width=205,height=270,units="mm",dpi=300,bg="white")

audit <- data.table(
  section=c("3.7","3.7","3.7","3.7","3.7","3.8","3.8","3.8","3.8","3.8"),
  metric=c("47_gene_ARI","pure_proliferation_full_ARI","paired_delta_median","paired_delta_CI_low","paired_delta_CI_high",
    "external_CMS_cramers_V","within_CMS_min_Hedges_g","within_CMS_max_Hedges_g","increment_adjusted_R2_all","partial_R2_all"),
  value=c(con[analysis=="47_gene_deletion_sensitivity",ARI],con[analysis=="E2F_G2M_full",ARI],
    bootsum[contrast=="ARI_47_minus_E2F_G2M_full",delta_median],bootsum[contrast=="ARI_47_minus_E2F_G2M_full",CI95_low],
    bootsum[contrast=="ARI_47_minus_E2F_G2M_full",CI95_high],cmsg$cramers_V,min(within$hedges_g),max(within$hedges_g),
    inc[analysis=="all",adjusted_R2_increment],inc[analysis=="all",partial_R2])
)
fwrite(audit,file.path(out,"Results_3.7_3.8_frozen_number_audit.tsv"),sep="\t")

writeLines(c(
  "Figure 5. Proliferation and CMS sensitivity analyses define the biological scope of the frozen state.",
  "(A) Concordance of the 47-gene deletion, pure-proliferation, and proliferation-residualized solutions with the frozen 68-gene K=2 labels.",
  "(B) Paired sample-bootstrap distributions for the ARI advantage of the 47-gene deletion solution over two E2F/G2M controls.",
  "(C) Mean TCGA silhouette widths for the frozen, deletion, and pure-proliferation feature spaces.",
  "(D) K=2-6 consensus stability after residualizing the 68-gene matrix on E2F/G2M scores. PAC is the proportion of ambiguous consensus values; lower values indicate greater stability.",
  "(E) K=2 composition across independently available external CMS labels (n=503).",
  "(F) Within-CMS C2-versus-C1 metabolic-intensity effect sizes among CMSclassifier-evaluable samples (n=645). Error bars are 2,000-bootstrap 95% CIs.",
  "(G) Adjusted R2 for the covariate model before and after adding the continuous state score, in all evaluable samples and after excluding 13 source/vial-associated cases."
),file.path(out,"Figure5_legend_draft.md"))

capture.output(sessionInfo(),file=file.path(out,"R_sessionInfo.txt"))
cat("Figure 5 completed\n")
