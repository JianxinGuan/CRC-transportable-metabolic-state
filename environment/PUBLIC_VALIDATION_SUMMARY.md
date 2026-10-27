# Numerical validation summary

Seven core chains were rerun in a clean ASCII-path test directory and compared with frozen outputs. The reproduced checks covered the Meta-CRISPR intersection, TCGA normalization and consensus clustering, fixed-centroid external projection, 624-patient acidification/metabolic models, CMS/MSI and nested CMS models, survival/immune outputs, and GSE132465 single-cell processing.

The validation did not rerun raw preprocessing for every GEO accession and did not establish a single-command end-to-end reconstruction of every figure.

Key reproduced checkpoints included:

- 1,647 robust and 1,500 strict human meta-analysis candidates;
- 174 strict one-to-one Meta-CRISPR candidates and the frozen 68-gene subset;
- TCGA K=2 group sizes of 214 and 433;
- GSE39582 273/293 and GSE17536 87/90 fixed-centroid group sizes;
- GSE39582 and GSE17536 mean silhouettes of 0.1034438652 and 0.0907423697;
- GSE132465: 63,689 cells, 23 patients, 68/68 genes, paired Wilcoxon P = 2.384185791015625e-06.

Exact software sessions are provided in the accompanying `sessionInfo_*.txt` files.
