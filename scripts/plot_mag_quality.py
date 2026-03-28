#!/usr/bin/env python3
"""
plot_mag_quality.py — Publication-quality figures for MAG analysis.

Outputs:
  - results/figures/mag_quality_barplot.pdf
  - results/figures/taxonomy_piechart.pdf
  - results/figures/abundance_heatmap.pdf

Usage:
    python3 scripts/plot_mag_quality.py \\
        --input results/summary_table.tsv \\
        --outdir results/figures
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

try:
    import pandas as pd
except ImportError:
    print("ERROR: pandas required. Install: pip install pandas matplotlib numpy")
    sys.exit(1)

QUALITY_COLORS = {
    "HQ": "#2ecc71",    # Green
    "MQ": "#3498db",    # Blue
    "LQ": "#f39c12",    # Orange
    "Poor": "#e74c3c",  # Red
    "Unknown": "#95a5a6" # Grey
}

def load_data(path):
    """Load summary table TSV."""
    if not os.path.exists(path):
        print(f"ERROR: Summary table not found: {path}", file=sys.stderr)
        sys.exit(1)
    df = pd.read_csv(path, sep="\t")
    # Clean mag_id
    df["mag_id"] = df["mag_id"].astype(str)
    return df

def plot_quality_barplot(df, outdir):
    """
    Bar chart: Completeness (blue) and Contamination (red) per MAG,
    colored by quality tier.
    """
    if df.empty:
        print("WARNING: No data for barplot", file=sys.stderr)
        return

    df = df.copy()
    df = df.sort_values("completeness", ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(max(10, len(df) * 0.5), 5))

    n = len(df)
    x = np.arange(n)
    bar_width = 0.6

    colors = [QUALITY_COLORS.get(q, "#95a5a6") for q in df["quality"]]

    bars = ax.bar(x, df["completeness"], bar_width,
                  color=colors, edgecolor="white", linewidth=0.5,
                  label="Completeness", alpha=0.85)

    # Contamination as a transparent overlay
    ax.bar(x, df["contamination"], bar_width,
           color="#e74c3c", alpha=0.4, label="Contamination")

    # Threshold lines
    ax.axhline(y=90, color="#2ecc71", linestyle="--", linewidth=1.5, alpha=0.8, label="HQ threshold (90%)")
    ax.axhline(y=50, color="#3498db", linestyle="--", linewidth=1.5, alpha=0.8, label="MQ/LQ threshold (50%)")
    ax.axhline(y=5, color="#e74c3c", linestyle=":", linewidth=1.5, alpha=0.8, label="Contam. threshold (5%)")

    ax.set_xlabel("MAG", fontsize=11)
    ax.set_ylabel("Percentage (%)", fontsize=11)
    ax.set_title("MAG Quality: Completeness & Contamination (MIMAG Standard)", fontsize=12, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(df["mag_id"], rotation=90, fontsize=7)
    ax.set_ylim(0, 105)
    ax.legend(loc="upper right", fontsize=8)

    # Legend patches for quality tiers
    handles = [
        mpatches.Patch(color=QUALITY_COLORS["HQ"], label="HQ (≥90% comp, <5% contam)"),
        mpatches.Patch(color=QUALITY_COLORS["MQ"], label="MQ (≥50% comp, <10% contam)"),
        mpatches.Patch(color=QUALITY_COLORS["LQ"], label="LQ (≥50% comp, <5% contam)"),
        mpatches.Patch(color=QUALITY_COLORS["Poor"], label="Poor (<50% comp or ≥10% contam)"),
    ]
    ax.legend(handles=handles + ax.get_legend_handles_labels()[0], loc="upper right", fontsize=7)

    plt.tight_layout()
    out = os.path.join(outdir, "mag_quality_barplot.pdf")
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  ✓ Saved: {out}")

def plot_taxonomy_piechart(df, outdir):
    """Phylum-level pie chart of MAGs."""
    if df.empty or "phylum" not in df.columns:
        return

    phylum_counts = df["phylum"].fillna("Unknown").value_counts()

    # Shorten long labels
    labels = [p.replace("p__", "")[:30] for p in phylum_counts.index]

    fig, ax = plt.subplots(figsize=(8, 8))
    colors = plt.cm.Set3(np.linspace(0, 1, len(phylum_counts)))

    wedges, texts, autotexts = ax.pie(
        phylum_counts,
        labels=labels,
        autopct=lambda p: f"{p:.1f}%\n({int(p/100*phylum_counts.sum())})" if p > 3 else "",
        colors=colors,
        startangle=90,
        pctdistance=0.75,
        textprops={"fontsize": 8},
    )
    for t in texts:
        t.set_fontsize(7)

    ax.set_title("MAG Taxonomy Distribution (Phylum Level)", fontsize=12, fontweight="bold")
    plt.tight_layout()
    out = os.path.join(outdir, "taxonomy_piechart.pdf")
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  ✓ Saved: {out}")

def plot_abundance_heatmap(df, outdir):
    """Sample × MAG abundance heatmap."""
    if df.empty or "abundance" not in df.columns:
        print("WARNING: No abundance data for heatmap", file=sys.stderr)
        return

    df = df.copy()
    df = df.sort_values("abundance", ascending=False).reset_index(drop=True)

    # Select top 30 MAGs for readability
    if len(df) > 30:
        df = df.head(30)

    fig, ax = plt.subplots(figsize=(10, max(6, len(df) * 0.3)))

    abundances = df[["mag_id", "abundance", "genus"]].copy()
    abundances = abundances.set_index("mag_id")

    matrix = np.log10(abundances["abundance"].replace(0, np.nan).fillna(1e-10) + 1).values.reshape(-1, 1)

    im = ax.imshow(matrix, aspect="auto", cmap="YlOrRd", interpolation="nearest")

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels([f"{row['mag_id']} ({row['genus']})" for _, row in df.iterrows()], fontsize=7)
    ax.set_xticks([0])
    ax.set_xticklabels(["Abundance (log10)"], fontsize=9)
    ax.set_title("MAG Abundance (Top 30 by Coverage)", fontsize=12, fontweight="bold")

    plt.colorbar(im, ax=ax, label="log10(abundance + 1)", shrink=0.5)
    plt.tight_layout()
    out = os.path.join(outdir, "abundance_heatmap.pdf")
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  ✓ Saved: {out}")

def main():
    parser = argparse.ArgumentParser(description="Generate MAG analysis figures")
    parser.add_argument("--input", required=True, help="Path to summary_table.tsv")
    parser.add_argument("--outdir", default="results/figures", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = load_data(args.input)

    print(f"\nGenerating figures from {args.input} ({len(df)} MAGs)...")

    plot_quality_barplot(df, args.outdir)
    plot_taxonomy_piechart(df, args.outdir)
    plot_abundance_heatmap(df, args.outdir)

    print(f"\n✓ All figures saved to: {args.outdir}/")

if __name__ == "__main__":
    main()
