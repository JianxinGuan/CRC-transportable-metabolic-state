# Public-release candidate status

Prepared on 2026-10-27 from the private provenance repository without copying its Git history.

## Completed P0 cleanup

- The supported run path contains only ten numerically tested core scripts.
- Fourteen unvalidated historical/figure scripts are isolated under `archive/provenance_only_scripts/`.
- The obsolete internal freeze script was excluded.
- Internal path/secret scan reports and private release checklists were excluded.
- README, script inventory, run order, expected outputs, and data-access metadata were rewritten in UTF-8.
- Seven validated chains are stated explicitly without claiming raw-data end-to-end reproduction.
- Source-data URLs and redistribution boundaries are documented.
- Python dependencies are pinned in `requirements.txt`.
- Supporting Figures S1-S3 and `S1_Table.xlsx` are byte-identical to the current journal-upload directory.
- Both PDF and PNG versions of Supporting Figures S1-S3 are intentionally retained: PDF for vector review/production and PNG for direct preview.
- Core Python scripts compile successfully; core R scripts parse successfully when copied to an ASCII-only test path.
- No manuscript, raw matrix, credential, token, archive, or restricted-source object is included.
- The manuscript Methods citation for the GSE132465 summary was corrected from sheets `S12a-S12d` to sheet `S13` in `CRC_three_phase_workflow_PLOS_ONE_contact_code_S13_audited.docx`.

## Remaining blockers before public upload

1. Create and review the new clean-history Git repository from this directory. Do not reuse the private provenance repository history without an explicit replacement decision.
2. Decide whether the clean history will be pushed to a new public repository or will replace the existing private remote history.
