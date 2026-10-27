import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans


def main():
    parser = argparse.ArgumentParser(description="Reproduce stage 24 custom consensus k-means.")
    parser.add_argument("--zscore", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--seed", type=int, default=20260827)
    parser.add_argument("--runs", type=int, default=500)
    parser.add_argument("--sample-fraction", type=float, default=0.8)
    parser.add_argument("--gene-fraction", type=float, default=0.8)
    parser.add_argument("--k-min", type=int, default=2)
    parser.add_argument("--k-max", type=int, default=6)
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.zscore, sep="\t", index_col=0)
    matrix = frame.to_numpy(float).T
    samples = frame.columns.tolist()
    if frame.shape != (68, 647):
        raise ValueError(f"Expected 68 genes x 647 samples, observed {frame.shape}")

    rng = np.random.default_rng(args.seed)
    summaries = []
    for k in range(args.k_min, args.k_max + 1):
        co_cluster = np.zeros((len(samples), len(samples)), dtype=float)
        co_sample = np.zeros_like(co_cluster)
        for bootstrap_index in range(args.runs):
            sample_index = rng.choice(
                len(samples), size=max(2, int(args.sample_fraction * len(samples))), replace=False
            )
            gene_index = rng.choice(
                matrix.shape[1], size=max(2, int(args.gene_fraction * matrix.shape[1])), replace=False
            )
            labels = KMeans(
                n_clusters=k, n_init=10, random_state=bootstrap_index
            ).fit_predict(matrix[sample_index][:, gene_index])
            for cluster in range(k):
                members = sample_index[labels == cluster]
                co_cluster[np.ix_(members, members)] += 1
            co_sample[np.ix_(sample_index, sample_index)] += 1

        consensus = np.divide(
            co_cluster, co_sample, out=np.zeros_like(co_cluster), where=co_sample > 0
        )
        np.fill_diagonal(consensus, 1)
        pd.DataFrame(consensus, index=samples, columns=samples).to_csv(
            out / f"consensus_matrix_k{k}.tsv.gz", compression="gzip", sep="\t"
        )

        final_fit = KMeans(n_clusters=k, n_init=50, random_state=1).fit(matrix)
        pd.DataFrame({"sample_id": samples, "cluster": final_fit.labels_ + 1}).to_csv(
            out / f"cluster_assignments_k{k}.tsv", sep="\t", index=False
        )
        upper = consensus[np.triu_indices_from(consensus, 1)]
        summaries.append({
            "k": k,
            "mean_consensus": float(upper.mean()),
            "sd_consensus": float(upper.std()),
            "min_consensus": float(upper.min()),
            "max_consensus": float(upper.max()),
            "cluster_sizes": ";".join(map(str, np.bincount(final_fit.labels_))),
        })

    pd.DataFrame(summaries).to_csv(out / "consensus_clustering_stability.tsv", sep="\t", index=False)
    best = max(summaries, key=lambda row: row["mean_consensus"] - row["sd_consensus"])
    (out / "recommended_k.txt").write_text(f"{best['k']}\n", encoding="ascii")
    (out / "README.txt").write_text(
        "500 resampling runs; 80% samples and 80% genes per run; KMeans n_init=10; "
        "final labels fit all samples with n_init=50. Recommended K maximizes "
        "mean_consensus - sd_consensus; inspect matrices and cluster sizes before "
        "biological interpretation.\n",
        encoding="ascii",
    )


if __name__ == "__main__":
    main()
