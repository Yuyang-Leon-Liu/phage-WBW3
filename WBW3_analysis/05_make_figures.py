#!/usr/bin/env python3
"""Step 5 - Redraw the data-driven figures from the result tables.

Only the panels that come straight from the CSV tables are redrawn here
(the phenotypic atlas, the marine-homolog histograms, and the CRISPR library
summary). The genome map, the phylogenetic trees and the protein-structure
figure need the upstream tools from steps 1-3 and are not redrawn by this script.

Run from the folder that contains the data/ directory:
    python 05_make_figures.py
Figures are written to the current folder as PNG.
"""
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["Liberation Sans", "DejaVu Sans"]

BLUE, ORANGE, GREEN, GREY = "steelblue", "darkorange", "seagreen", "gray"


# English labels for the five sampling sites (the raw table uses Chinese names)
SAMPLE_EN = {
    "鸡场污水": "chicken-farm sewage",
    "菜市场污水": "market sewage",
    "鸭子粪便": "duck feces",
    "猪场污水": "pig-farm sewage",
    "河水": "river water",
}


def phenotypic_atlas():
    """Figure 2 - the 54-phage collection at a glance."""
    iso = pd.read_csv("data/isolation_stats.csv")
    iso["sample"] = iso["sample"].map(lambda s: SAMPLE_EN.get(s, s))
    atlas = pd.read_csv("data/phage_phenotypic_atlas.csv")

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))

    # A - how many distinct phages per sample, and the isolation rate
    ax = axes[0, 0]
    ax.bar(iso["sample"], iso["distinct_phages"], color=BLUE)
    ax.set_ylabel("distinct phages")
    ax.set_title("A  Isolation per sample", loc="left", fontweight="bold")
    ax.tick_params(axis="x", rotation=45)

    # B - lysate titer by lysis group
    ax = axes[0, 1]
    atlas.boxplot(column="titer_pfu_ml", by="lysis_group", ax=ax)
    fig.suptitle("")   # drop pandas' default figure title
    ax.set_title("")   # drop pandas' column-name title
    ax.set_xlabel("")  # drop the "lysis_group" axis label
    ax.grid(False)
    ax.set_yscale("log")
    ax.set_ylabel("titer (PFU/mL)")
    ax.set_title("B  Lysate titer by group", loc="left", fontweight="bold")
    ax.tick_params(axis="x", rotation=45)

    # C - virion size by morphotype
    ax = axes[1, 0]
    for morph, grp in atlas.groupby("morphotype"):
        ax.scatter(grp["head_width_nm_v"], grp["tail_length_nm_v"], label=morph, s=18)
    ax.set_xlabel("head width (nm)")
    ax.set_ylabel("tail length (nm)")
    ax.set_title("C  Virion morphometrics", loc="left", fontweight="bold")
    ax.legend(fontsize=6)

    # D - how many phages in each lysis group
    ax = axes[1, 1]
    counts = atlas["lysis_group"].value_counts().sort_index()
    ax.bar(range(len(counts)), counts.values, color=GREEN)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("number of phages")
    ax.set_title("D  Six lysis-kinetic groups", loc="left", fontweight="bold")

    fig.tight_layout()
    fig.savefig("fig2_phenotypic_atlas.png", dpi=200)
    print("wrote fig2_phenotypic_atlas.png")


def marine_landscape():
    """Figure 6 (panels A and C) - the marine endolysin homologs."""
    hits = pd.read_csv("data/wbw3_vs_omd_all_hits.csv")
    lys = pd.read_csv("data/marine_endolysin_candidates.csv")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # A - identity distribution of the marine endolysin homologs
    ax = axes[0]
    ax.hist(lys["pident_pct"], bins=20, color=BLUE, edgecolor="white")
    ax.axvline(50, color=ORANGE, linestyle="--", label="50% identity")
    ax.set_xlabel("amino-acid identity to WBW3 endolysin (%)")
    ax.set_ylabel("marine homologs")
    ax.set_title("A  235 divergent marine endolysins", loc="left", fontweight="bold")
    ax.legend()

    # C - marine homologs per functional category
    ax = axes[1]
    counts = hits["function"].value_counts()
    colors = [ORANGE if c == "lysis" else GREY for c in counts.index]
    ax.barh(range(len(counts)), counts.values, color=colors)
    ax.set_yticks(range(len(counts)))
    ax.set_yticklabels(counts.index, fontsize=7)
    ax.set_xlabel("marine homologs")
    ax.set_title("C  Marine homologs by category", loc="left", fontweight="bold")

    fig.tight_layout()
    fig.savefig("fig6_marine_panels.png", dpi=200)
    print("wrote fig6_marine_panels.png")


def crispr_library():
    """Figure 6 (panels D and E) - the essential-gene guide library."""
    g = pd.read_csv("data/essential_gene_gRNA_categorized.csv")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # D - essential genes per category
    ax = axes[0]
    counts = g["category"].value_counts()
    ax.barh(range(len(counts)), counts.values, color=BLUE)
    ax.set_yticks(range(len(counts)))
    ax.set_yticklabels(counts.index, fontsize=8)
    ax.set_xlabel("essential genes")
    ax.set_title("D  Essential-gene target landscape", loc="left", fontweight="bold")

    # E - dual-strand design (one guide per strand; 320 of 324 genes got both)
    ax = axes[1]
    n_top = g["gRNA_top"].notna().sum()
    n_bot = g["gRNA_bottom"].notna().sum()
    ax.pie([n_top, n_bot], labels=["top strand", "bottom strand"],
           colors=[BLUE, ORANGE], autopct="%d", startangle=90,
           wedgeprops=dict(width=0.4))
    ax.set_title("E  Dual-strand gRNA design", loc="left", fontweight="bold")

    fig.tight_layout()
    fig.savefig("fig6_crispr_panels.png", dpi=200)
    print("wrote fig6_crispr_panels.png")


if __name__ == "__main__":
    phenotypic_atlas()
    marine_landscape()
    crispr_library()
