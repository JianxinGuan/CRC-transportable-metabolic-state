import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Reproduce stage 23 CRC 68-gene normalization.")
    parser.add_argument("--genes", required=True, help="One gene symbol per line, in frozen order.")
    parser.add_argument("--tpm", required=True, help="TCGA primary-tumor TPM TSV or TSV.GZ.")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(args.genes, header=None, sep="\t")[0].astype(str).tolist()
    if len(genes) != 68 or len(set(g.upper() for g in genes)) != 68:
        raise ValueError("The frozen gene list must contain 68 unique symbols.")

    expression = pd.read_csv(args.tpm, sep="\t", compression="infer")
    required = {"gene_id", "gene_name", "gene_type"}
    if not required.issubset(expression.columns):
        raise ValueError(f"Missing required columns: {sorted(required - set(expression.columns))}")
    expression["_symbol_key"] = expression["gene_name"].astype(str).str.upper()
    expression = expression.drop_duplicates("_symbol_key", keep="first").set_index("_symbol_key")
    selected = expression.reindex([g.upper() for g in genes])
    missing = [genes[i] for i, value in enumerate(selected["gene_name"].isna()) if value]
    if missing:
        raise ValueError(f"Frozen genes absent from TPM matrix: {missing}")

    values = selected.iloc[:, 3:].apply(pd.to_numeric, errors="coerce").fillna(0)
    log2_tpm = np.log2(values + 1)
    row_mean = log2_tpm.mean(axis=1)
    row_sd = log2_tpm.std(axis=1, ddof=1)
    zscore = log2_tpm.sub(row_mean, axis=0).div(row_sd.replace(0, np.nan), axis=0).fillna(0)
    log2_tpm.index = genes
    zscore.index = genes

    log2_tpm.to_csv(out / "CRC68_TCGA_log2TPM.tsv", sep="\t")
    zscore.to_csv(out / "CRC68_TCGA_gene_zscore.tsv", sep="\t")
    pd.DataFrame({
        "gene_symbol": genes,
        "mean_log2TPM": log2_tpm.mean(axis=1).values,
        "sd_log2TPM": log2_tpm.std(axis=1, ddof=1).values,
    }).to_csv(out / "CRC68_normalization_audit.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
