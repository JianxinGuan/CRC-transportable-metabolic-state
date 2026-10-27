invisible(Sys.setlocale("LC_CTYPE", "Chinese"))
suppressPackageStartupMessages({library(data.table); library(clusterProfiler); library(org.Hs.eg.db); library(igraph)})
out <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "geo_stage17")
dir.create(file.path(out,"01_functional_enrichment"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,"02_PPI"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,"03_CRC_expression_filter"),recursive=TRUE,showWarnings=FALSE)

# Input: strict Meta candidates intersected with Ensembl 113 one-to-one orthologs.
dual <- fread(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "geo_stage16/04_crispr_intersections/prioritized_strict_one_to_one_dual_evidence_candidates.tsv"))
genes <- unique(toupper(dual$gene_symbol))
stopifnot(length(genes)==174L)
fwrite(data.table(gene_symbol=genes),file.path(out,"dual_evidence_180_genes.tsv"),sep="\t")

# GO enrichment (all 180 genes; BP/MF/CC) with human genome background.
eg <- mapIds(org.Hs.eg.db, keys=genes, column="ENTREZID", keytype="SYMBOL", multiVals="first")
eg <- unique(na.omit(unname(eg)))
go <- enrichGO(gene=eg, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="ALL",
               pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1, readable=TRUE,
               universe=keys(org.Hs.eg.db, keytype="ENTREZID"))
go_df <- as.data.table(as.data.frame(go))
if(nrow(go_df)) fwrite(go_df,file.path(out,"01_functional_enrichment","GO_ALL_180.tsv"),sep="\t")
go_sig <- go_df[p.adjust < 0.05]
if(nrow(go_sig)) fwrite(go_sig,file.path(out,"01_functional_enrichment","GO_ALL_180_FDR05.tsv"),sep="\t")
writeLines(c(paste0("input_genes=",length(genes)),paste0("mapped_entrez=",length(eg)),
             paste0("GO_terms_tested=",nrow(go_df)),paste0("GO_terms_FDR05=",nrow(go_sig))),
           file.path(out,"01_functional_enrichment","enrichment_summary.txt"))

# STRING v12 network, score >=700, using the public API and release metadata.
query <- paste(genes, collapse="%0d")
url <- paste0("https://string-db.org/api/tsv/network?identifiers=",query,
              "&species=9606&required_score=700&caller_identity=CRC_acidosis_stage17")
tf <- tempfile(fileext=".tsv")
status <- tryCatch(download.file(url,tf,mode="wb",quiet=TRUE),error=function(e) 1L)
if(status==0L && file.exists(tf) && file.info(tf)$size>0) {
  ppi <- fread(tf,fill=TRUE)
  fwrite(ppi,file.path(out,"02_PPI","STRING_network_score700.tsv"),sep="\t")
  if(nrow(ppi)) {
    edges <- ppi[,.(from=preferredName_A,to=preferredName_B,score=score)]
    g <- graph_from_data_frame(edges,directed=FALSE)
    deg <- sort(degree(g),decreasing=TRUE)
    hubs <- data.table(gene_symbol=names(deg),degree=as.integer(deg))
    fwrite(hubs,file.path(out,"02_PPI","STRING_degree_hubs.tsv"),sep="\t")
    fwrite(data.table(nodes=vcount(g),edges=ecount(g),connected_components=components(g)$no),
           file.path(out,"02_PPI","PPI_summary.tsv"),sep="\t")
  }
} else writeLines("STRING API unavailable; no PPI fabricated.",file.path(out,"02_PPI","PPI_status.txt"))

# CRC expression audit. These RDS files are data.frames with rows=samples,
# columns=gene symbols plus OS.status/OS.time; no sample-level TCGA matrix is
# currently present, so TCGA is explicitly audited as pending.
rds <- list.files(file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "crc_external"),pattern="_expr_surv\\.rds$",full.names=TRUE)
keep <- setdiff(genes,c("OS.STATUS","OS.TIME"))
rows <- list()
for(f in rds){
  z <- readRDS(f); n <- nrow(z)
  present <- intersect(keep,colnames(z));
  if(length(present)) {
    mat <- as.matrix(z[,present,drop=FALSE]); storage.mode(mat)<-"double"
    detect <- colMeans(is.finite(mat) & mat != 0)
    vars <- apply(mat,2,var,na.rm=TRUE)
    rows[[length(rows)+1L]] <- data.table(dataset=sub("_expr_surv.rds$","",basename(f)),
      n_samples=n, n_genes_present=length(present),
      gene_symbol=present, detection_rate=as.numeric(detect), variance=as.numeric(vars))
  }
  rm(z); gc()
}
expr_audit <- rbindlist(rows,fill=TRUE)
expr_audit[, variance_percentile_within_180 := frank(variance, ties.method="average")/.N,
           by=dataset]
fwrite(expr_audit,file.path(out,"03_CRC_expression_filter","CRC_external_gene_detection_variance.tsv.gz"),sep="\t",compress="gzip")

# Predefined, unsupervised filter: present in >=5/8 external cohorts, detection
# >=80% within each contributing cohort, and positive variance in >=5 cohorts.
cohort_stats <- expr_audit[,.(cohorts_present=sum(detection_rate>=0.8),
                              median_detection=median(detection_rate[detection_rate>=0.8],na.rm=TRUE),
                              cohorts_variable=sum(is.finite(variance)&variance>1e-8),
                              high_variance_cohorts=sum(variance_percentile_within_180>=0.5,na.rm=TRUE),
                              median_variance_percentile=median(variance_percentile_within_180,na.rm=TRUE),
                              median_variance=median(variance[is.finite(variance)],na.rm=TRUE)),by=gene_symbol]
cohort_stats[, external_crc_pass := cohorts_present>=5 & cohorts_variable>=5]
fwrite(cohort_stats,file.path(out,"03_CRC_expression_filter","CRC_external_gene_filter_summary.tsv"),sep="\t")

tcga_dir <- file.path(Sys.getenv("CRC_ARCHIVE_ROOT", unset=""), "D", "tcga_coadread")
tcga_expr <- list.files(tcga_dir,recursive=TRUE,full.names=TRUE,pattern="(expression|counts|FPKM|TPM).*(tsv|csv|rds|gz)$")
tcga_status <- data.table(tcga_expression_files=length(tcga_expr),tcga_expression_available=length(tcga_expr)>0,
                           note=if(length(tcga_expr)==0) "TCGA clinical/MSI/CMS present; expression matrix must be added before final TCGA expression filter." else paste(tcga_expr,collapse=";"))
fwrite(tcga_status,file.path(out,"03_CRC_expression_filter","TCGA_expression_availability.tsv"),sep="\t")

# Functional refinement is deliberately independent of clustering. Genes are
# retained if they pass external CRC expression criteria; all exclusions are logged.
final <- data.table(gene_symbol=genes)
final <- merge(final,cohort_stats,by="gene_symbol",all.x=TRUE)
final <- merge(final,dual[,.(gene_symbol=toupper(gene_symbol),mouse_symbol,support_n,
  random_effect,random_FDR,k,direction_consistency,priority_tier)],by="gene_symbol",all.x=TRUE)
final[, strong_dual_evidence := priority_tier %chin% c("Tier_1","Tier_2")]
final[, reproducibly_variable := high_variance_cohorts>=4]
final[, final_crc_subtype_input := external_crc_pass == TRUE & strong_dual_evidence & reproducibly_variable]
final[, exclusion_reason := fcase(
  final_crc_subtype_input, "retained_provisional_CRC_subtype_input",
  !external_crc_pass, "excluded_external_CRC_detection_or_variance",
  !strong_dual_evidence, "excluded_Tier3_dual_evidence",
  !reproducibly_variable, "excluded_not_top50pct_variable_in_4_cohorts",
  default="excluded_other")]
fwrite(final,file.path(out,"03_CRC_expression_filter","CRC_contextualized_180_gene_audit.tsv"),sep="\t")
fwrite(final[final_crc_subtype_input==TRUE,.(gene_symbol)],
       file.path(out,"CRC_final_subtype_gene_set_external_CRC_only.txt"),col.names=FALSE)

summary <- data.table(stage=c("dual_evidence_input","external_CRC_availability_pass",
 "strong_Tier1_or_Tier2","reproducibly_variable","provisional_subtype_gene_set"),
 n=c(nrow(final),sum(final$external_crc_pass,na.rm=TRUE),sum(final$strong_dual_evidence,na.rm=TRUE),
     sum(final$reproducibly_variable,na.rm=TRUE),sum(final$final_crc_subtype_input,na.rm=TRUE)),
 rule=c("180 strict Meta x CRISPR one-to-one","present >=5 cohorts, detection >=80%, variable >=5 cohorts",
        "Stage16 Tier1 or Tier2","top 50% variance among 180 in >=4 of 8 CRC cohorts",
        "all three filters; provisional until TCGA expression verification"))
fwrite(summary,file.path(out,"03_CRC_expression_filter","CRC_filter_summary.tsv"),sep="\t")
writeLines(c(
 "This is CRC-contextualized gene-set preparation before clustering.",
 "The external CRC filter is prespecified and independent of subtype labels.",
 "Provisional subtype input additionally requires Stage16 Tier1/2 and top-50% variance in at least 4 of 8 CRC cohorts.",
 "TCGA-COAD/READ expression matrix was not found in the designated TCGA directory; final TCGA expression filtering remains pending.",
 "GSE cell line and mixed tissue cohorts are not used here; only CRC survival RDS cohorts are used for expression availability.",
 "GO and STRING results are prioritization evidence, not causal validation.",
 "PPI source: STRING public API, species 9606, required_score=700; record API response and failure status."
),file.path(out,"README_STAGE17.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"R_sessionInfo.txt"))
cat("Stage17 completed; external CRC retained:",sum(final$final_crc_subtype_input,na.rm=TRUE),"of",nrow(final),"\n")
