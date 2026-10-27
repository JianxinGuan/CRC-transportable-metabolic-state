import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


CELL_ORDER = [
    "Epithelial cells",
    "Stromal cells",
    "T cells",
    "B cells",
    "Myeloids",
    "Mast cells",
]


def main():
    parser = argparse.ArgumentParser(description="Generate Supporting Figure S1 for GSE132465.")
    parser.add_argument("--cell-scores", required=True)
    parser.add_argument("--paired-scores", required=True)
    parser.add_argument("--test-results", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    cells = pd.read_csv(args.cell_scores, sep="\t")
    paired = pd.read_csv(args.paired_scores, sep="\t")
    test = pd.read_csv(args.test_results, sep="\t").iloc[0]
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    missing_types = set(CELL_ORDER) - set(cells["Cell_type"].dropna().unique())
    if missing_types:
        raise ValueError(f"Missing expected cell types: {sorted(missing_types)}")
    if int(test["n_patients_paired"]) != len(paired):
        raise ValueError("Paired-patient count does not match the frozen test table")

    # Deterministic visual subsampling only; boxplots use all cells.
    point_data = pd.concat(
        [
            group.sample(min(len(group), 2000), random_state=132465)
            for _, group in cells.groupby("Cell_type", sort=False)
        ],
        ignore_index=True,
    )

    sns.set_theme(style="whitegrid", context="talk")
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={"width_ratios": [1.25, 1]})

    sns.boxplot(
        data=cells,
        x="Cell_type",
        y="score68",
        order=CELL_ORDER,
        color="#7FA6BF",
        width=0.65,
        showfliers=False,
        linewidth=1.5,
        ax=axes[0],
    )
    sns.stripplot(
        data=point_data,
        x="Cell_type",
        y="score68",
        order=CELL_ORDER,
        color="#456A82",
        alpha=0.16,
        size=1.5,
        jitter=0.28,
        ax=axes[0],
    )
    axes[0].set_title("Cell-level score distribution")
    axes[0].set_xlabel("Cell type")
    axes[0].set_ylabel("68-gene score")
    axes[0].tick_params(axis="x", rotation=38)

    x_left, x_right = 0, 1
    for _, row in paired.iterrows():
        axes[1].plot(
            [x_left, x_right],
            [row["malignant"], row["non_malignant"]],
            color="#B9C4CF",
            linewidth=1.1,
            alpha=0.75,
            zorder=1,
        )
    axes[1].scatter(np.repeat(x_left, len(paired)), paired["malignant"], color="#C84F4F", s=42, zorder=2)
    axes[1].scatter(np.repeat(x_right, len(paired)), paired["non_malignant"], color="#4B78B8", s=42, zorder=2)
    axes[1].set_xticks([0, 1], ["Malignant/\nepithelial", "Non-malignant"])
    axes[1].set_xlim(-0.35, 1.35)
    axes[1].set_title("Paired patient-level comparison")
    axes[1].set_ylabel("Patient-level median 68-gene score")
    p_value = float(test["P_value"])
    exponent = int(np.floor(np.log10(p_value)))
    coefficient = p_value / (10 ** exponent)
    annotation = rf"n = {len(paired)}; Wilcoxon $P = {coefficient:.2f} \times 10^{{{exponent}}}$"
    axes[1].text(0.5, 0.97, annotation, ha="center", va="top", transform=axes[1].transAxes, fontsize=14)

    fig.suptitle("GSE132465 primary CRC single-cell validation", fontsize=22, fontweight="bold", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.96), w_pad=2.2)
    fig.savefig(output_dir / "S1_Fig.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(output_dir / "S1_Fig.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
