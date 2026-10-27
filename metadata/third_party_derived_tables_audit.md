# Third-party-derived summary-table audit

Audit date: 2026-10-27. This review concerns the eight derived TSV files included in `results/summary_tables/`. No source matrices, PDB coordinate files, GDSC training objects, or database dumps are included.

| File | Source dependency | What is included | Reviewer value | Public-release decision | Required attribution/condition |
|---|---|---|---|---|---|
| `S11a_drug_screening_funnel.tsv` | GDSC2/oncoPredict | Candidate counts and filtering stages only | Documents the exploratory funnel | Keep as supplementary, non-core | Cite GDSC2 and oncoPredict; state that training objects are not redistributed |
| `S11b_final8_predicted_sensitivity.tsv` | GDSC2/oncoPredict plus project-derived scores | Eight predicted-sensitivity summaries and derived ranking fields | Allows reviewers to inspect exploratory numbers | Keep as supplementary, non-core | Cite source resources; label predictions as hypothesis-generating and not measured patient sensitivity |
| `S12a_direct_binding_and_cocrystal_audit.tsv` | RCSB PDB and literature | PDB identifiers, target/ligand interpretation, resolution, and audit decisions | Supports structural interpretation boundaries | Keep as supplementary, non-core | Cite PDB entries 5V5Y and 9D0P; no coordinates are redistributed |
| `S12b_WEE1_PLK1_structure_audit.tsv` | RCSB PDB | Structure metadata and SHA-256 values for source coordinate files | Makes the docking provenance inspectable | Keep as supplementary, non-core | Cite RCSB PDB; hashes identify source files and do not replace the source license |
| `S12c_redocking_RMSD.tsv` | PDB structures and docking output | Docking poses, energies, RMSD, and pass/fail labels | Reports structural evidence and failure boundary | Keep as supplementary, non-core | Cite PDB and docking software; do not interpret as direct binding or efficacy |
| `S12d_interaction_retention.tsv` | PDB structures and PLIP-derived contacts | Contact-retention summaries and Jaccard values | Documents pose-quality comparison | Keep as supplementary, non-core | Cite PDB and PLIP; retain exploratory wording |
| `S14c_STRING_network_score400.tsv` | STRING database | Protein-pair identifiers, names, scores, and evidence-channel fields | Supports functional-network context for the 68-gene feature | Keep if the final manuscript uses STRING results | Cite STRING; retain the database version/access date and do not imply ownership of the database |
| `S4b_Ensembl113_complete_mapping_audit.tsv` | Ensembl release 113 and project-generated orthology mapping | Mouse-human mapping, support counts, LFC/FDR summaries, and inclusion flags | Essential audit trail for cross-species feature construction | Keep | Cite Ensembl release 113; state that this is a derived mapping audit, not a copy of the Ensembl database |

## Overall decision

None of these files should be treated as primary evidence. S11 and S12 are exploratory and can be moved to a supplementary-only directory if repository size or licensing review requires it. S4b is methodologically important and should remain in the public release. S14c should remain only if the final manuscript retains the STRING network analysis.

## Conditions before public release

1. Retain S11a/S11b as study-generated exploratory summaries. Do not add the GDSC training matrices, complete fitted-response files, oncoPredict internal objects, or the complete drug-by-sample prediction matrix.
2. Add source version/access date and citations for STRING and Ensembl to `metadata/data_sources.tsv`.
3. Keep PDB identifiers and derived structure metrics, but do not add PDB coordinate files or third-party database dumps.
4. In README and the manuscript, describe S11/S12 as exploratory and do not claim direct binding, clinical efficacy, or patient sensitivity.
5. Recompute the final manifest and SHA-256 file after any inclusion/exclusion decision.

## GDSC/oncoPredict final determination

The official GDSC FTP exposes GDSC2 raw and fitted-response release files as public downloads. The CancerRxGene site did not provide a clearly retrievable license statement specifically addressing redistribution of study-generated derived summaries. CRAN identifies oncoPredict version 1.3.1 as GPL-2. S11a and S11b contain only study-generated aggregate counts and an eight-candidate prediction summary; they do not contain source training matrices, complete fitted-response files, package code, or a complete prediction matrix. They are therefore retained with attribution and explicit exploratory/non-clinical wording. This practical publication assessment is not a transfer of ownership or a legal opinion.
