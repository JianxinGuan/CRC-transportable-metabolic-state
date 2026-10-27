# CRC transportable metabolic-proliferative state

Code and legally redistributable derived materials for a three-phase colorectal-cancer workflow: feature construction, fixed-centroid transportability assessment, and explicit boundary setting.

## Repository scope

This release contains frozen gene sets, ten numerically tested core scripts, derived summary tables, Supporting Information, software-session records, and integrity manifests. It excludes source expression matrices, controlled or license-restricted objects, manuscript files, editorial records, credentials, and local working archives.

The reported state is a continuous metabolic-proliferative expression state with lactate-related transcriptional support. The analyses do not establish direct intratumoral pH, acidosis specificity, causal adaptation, prognostic classification, or therapeutic efficacy.

## Reproducibility status

Seven core result chains were compared numerically with frozen analysis outputs:

1. Meta-CRISPR intersection and the 174-candidate evidence set;
2. TCGA 68-gene normalization and consensus clustering;
3. fixed-centroid projection into GSE39582 and GSE17536;
4. 624-patient acidification, lactate, pathway, and metabolic models;
5. CMS/MSI, within-CMS, and nested CMS models;
6. patient-level survival and immune sensitivity analyses;
7. GSE132465 single-cell scoring and patient-level comparison.

These tests used legally redistributable or locally staged frozen intermediate inputs. Raw preprocessing for every GEO accession and the complete figure-generation layer were not rerun as a single end-to-end public-layout workflow. Files in `archive/provenance_only_scripts/` are supplied only to document analytical provenance and are not part of the supported execution path.

## Contents

- `scripts/`: numerically tested core scripts;
- `archive/provenance_only_scripts/`: unvalidated historical/provenance scripts, excluded from the supported run path;
- `config/`: frozen parameters and gene sets;
- `metadata/`: data-access, script-status, run-order, and output-mapping documentation;
- `results/summary_tables/`: publication-facing derived summaries;
- `supporting_information/`: Supporting Figures and consolidated S1 Table workbook;
- `environment/`: software-session records and dependency guidance;
- `manifests/`: final file inventory and SHA-256 checksums.

## Source data

Source datasets are not redistributed. GEO accessions, GDC projects, and third-party resources are listed in `metadata/data_sources.tsv`; access and redistribution limits are described in `metadata/DATA_ACCESS.md`.

The GSE132465 script expects the following GEO processed files in a user-specified input directory:

- `GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz`
- `GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz`

## Running the supported scripts

Follow `metadata/run_order.md`. Each supported script requires explicit command-line paths or the documented staged input layout. Expected outputs and manuscript mappings are listed in `metadata/expected_outputs.tsv` and `metadata/script_output_manuscript_map.tsv`.

Example:

```bash
python scripts/05_single_cell/01_process_gse132465_primary_crc.py \
  --input-dir data/external/GSE132465 \
  --output-dir results/single_cell \
  --genes config/gene_sets/frozen_68_genes.tsv
```

## License and third-party materials

Repository source code is licensed under the MIT License. Supporting tables, figures, database-derived materials, and third-party content are not relicensed by the code license and remain subject to their source terms and the eventual article license.

## Citation

Citation metadata are provided in `CITATION.cff`. Cite the fixed release used for analysis rather than the moving default branch.
