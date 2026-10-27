import argparse
import gzip
import os

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

parser = argparse.ArgumentParser(description="Score the frozen 68-gene feature in GSE132465.")
parser.add_argument("--input-dir", required=True, help="Directory containing the two GEO processed .txt.gz files")
parser.add_argument("--output-dir", required=True, help="Output directory")
parser.add_argument("--genes", default="config/gene_sets/frozen_68_genes.tsv", help="TSV containing a gene column")
args = parser.parse_args()

base = os.path.abspath(args.input_dir)
out = os.path.abspath(args.output_dir)
os.makedirs(out, exist_ok=True)
annp = os.path.join(base, "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz")
matp = os.path.join(base, "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz")

gene_table = pd.read_csv(args.genes, sep="\t")
gene_col = "gene" if "gene" in gene_table.columns else gene_table.columns[0]
genes = gene_table[gene_col].astype(str).tolist()

ann = pd.read_csv(annp, sep="\t")
ann = ann.rename(columns={ann.columns[0]: "cell_id"})
with gzip.open(matp, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")[1:]
    rows = []
    for line in f:
        z = line.rstrip("\n").split("\t")
        if z[0] in genes:
            rows.append([z[0]] + z[1:])

expr = pd.DataFrame(rows, columns=["gene"] + header).set_index("gene").apply(pd.to_numeric)
common = [c for c in header if c in set(ann.cell_id)]
expr = expr[common]
ann = ann.set_index("cell_id").loc[common].reset_index()
lib = expr.sum(axis=0).replace(0, np.nan)
score = np.log1p(expr.div(lib, axis=1) * 1e4).mean(axis=0)
det = (expr > 0).sum(axis=0)
cell = ann.copy()
cell["score68"] = cell.cell_id.map(score)
cell["detected_68"] = cell.cell_id.map(det)
cell["detected_fraction"] = cell.detected_68 / len(genes)
cell.to_csv(os.path.join(out, "GSE132465_68gene_cell_scores.tsv"), sep="\t", index=False)

summary = cell.groupby(["Patient", "Sample", "Class", "Cell_type"], dropna=False).agg(
    n_cells=("cell_id", "size"),
    score_median=("score68", "median"),
    score_mean=("score68", "mean"),
    score_IQR=("score68", lambda x: x.quantile(.75) - x.quantile(.25)),
    detected_fraction_median=("detected_fraction", "median"),
).reset_index()
summary.to_csv(os.path.join(out, "GSE132465_68gene_patient_celltype_summary.tsv"), sep="\t", index=False)

pt = cell.groupby(["Patient", "Sample", "Class", "Cell_type"], dropna=False)["score68"].median().reset_index()
epithelial = pt[pt.Cell_type.str.contains("Epithelial|malignant", case=False, na=False)].groupby("Patient").score68.median()
other = pt[~pt.Cell_type.str.contains("Epithelial|malignant", case=False, na=False)].groupby("Patient").score68.median()
pair = pd.concat([epithelial.rename("epithelial_tumor_associated"), other.rename("other_cells")], axis=1).dropna()
pair["difference"] = pair.epithelial_tumor_associated - pair.other_cells
if len(pair) >= 3:
    w = wilcoxon(pair.epithelial_tumor_associated, pair.other_cells, alternative="two-sided")
    test = pd.DataFrame([{
        "n_patients_paired": len(pair),
        "median_epithelial_tumor_associated": pair.epithelial_tumor_associated.median(),
        "median_other_cells": pair.other_cells.median(),
        "median_difference": pair.difference.median(),
        "wilcoxon_W": w.statistic,
        "P_value": w.pvalue,
    }])
else:
    test = pd.DataFrame([{"n_patients_paired": len(pair), "P_value": np.nan}])
pair.to_csv(os.path.join(out, "GSE132465_patient_paired_epithelial_vs_other.tsv"), sep="\t")
test.to_csv(os.path.join(out, "GSE132465_patient_level_test.tsv"), sep="\t", index=False)

audit = pd.DataFrame([{
    "dataset": "GSE132465",
    "matrix_genes": expr.shape[0],
    "matrix_cells_original": len(header),
    "annotation_cells": len(ann),
    "matched_cells": len(common),
    "matched_fraction": len(common) / len(header),
    "n_patients": cell.Patient.nunique(),
    "n_samples": cell.Sample.nunique(),
    "n_cell_types": cell.Cell_type.nunique(),
    "genes_requested": len(genes),
    "genes_found": sum(g in expr.index for g in genes),
}])
audit.to_csv(os.path.join(out, "GSE132465_processing_summary.tsv"), sep="\t", index=False)
print(audit.to_string(index=False))
print(test.to_string(index=False))
