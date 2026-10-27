# Supported execution order

Only scripts listed as `release_ready` in `script_inventory.tsv` form the supported path. Public source data or locally staged frozen inputs are required and are not redistributed.

1. Run `scripts/01_feature_construction/04_meta_crispr_intersection.R` from frozen official-annotation Meta-ready summaries and the documented Ensembl 113 mapping input.
2. Run `scripts/02_tcga_state/01_normalize_68_gene_matrix.py`, followed by `02_consensus_clustering.py`.
3. Run `scripts/02_tcga_state/03_patient_dedup_survival.R` for the fixed 624-patient map and four survival endpoints.
4. Run `scripts/03_external_projection/01_gse39582_fixed_centroid_projection.R` for full-cohort fixed centroids and exclusion-13 sensitivity outputs.
5. Run the patient-level boundary scripts in this order:
   - `03_patient_cms_msi_metabolic.R`;
   - `04_exact_nested_models_624.R`;
   - `05_acidification_metabolic_624.R`;
   - `06_immune_sensitivity_624.R`.
6. Run `scripts/05_single_cell/01_process_gse132465_primary_crc.py` independently for GSE132465.

## Scope categories

- Primary: feature construction/freezing, full-cohort fixed-centroid projection, lactate-axis association, proliferation controls, and CMS increment model.
- Sensitivity: exclusion-13 centroids, patient deduplication and complete cases, alternate adjustment specifications, and proliferation controls.
- Exploratory: immune estimates, pharmacogenomic predictions, and structural analyses.

## Important limitation

The supported tests begin from frozen or locally staged intermediate inputs. The repository does not claim a single-command reconstruction from every raw GEO/GDC file.
