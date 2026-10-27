# Data access and redistribution

This repository does not redistribute source expression matrices, GDC source files, GDSC/oncoPredict training objects, or PDB coordinate files. Canonical access links are listed in `data_sources.tsv`.

## Public transcriptomic data

- Download GEO source or processed files from the accession pages listed in `data_sources.tsv`.
- Obtain TCGA-COAD and TCGA-READ files from the NCI Genomic Data Commons under the applicable GDC terms.
- The supported GSE132465 script expects:
  - `GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz`;
  - `GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz`.

## Third-party resources

- Obtain GDSC2/oncoPredict resources from their official sources under the applicable licenses.
- Obtain PDB structures 5V5Y and 9D0P directly from RCSB PDB.
- STRING and Ensembl-derived summary outputs retain their source attribution and are not presented as relicensed source databases.
- The RCSB PDB policy page states that wwPDB archive data files are available under the CC0 1.0 Public Domain Dedication; this repository nevertheless distributes only project-derived summaries, not coordinate files.
- S11a and S11b contain study-generated aggregate counts and an eight-candidate prediction summary. They do not redistribute GDSC source matrices, complete fitted-response files, oncoPredict internal objects, or a complete drug-by-sample prediction matrix. These summaries are exploratory and do not represent measured patient sensitivity.

## Included material

Only publication-facing derived summaries, frozen gene definitions, scripts, and Supporting Information intended for distribution are included. Users must retrieve source datasets independently. Inclusion of a derived table does not transfer ownership of or relicense its underlying database content.
