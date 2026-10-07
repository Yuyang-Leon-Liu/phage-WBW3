#!/usr/bin/env python3
"""Step 6 - Summarise the essential-gene CRISPR guide library.

The library holds one 23-nt NGG-PAM guide on each strand for 324 essential
E. coli genes. This script just tallies the library by functional category and
checks that every gene really has both a top- and a bottom-strand guide.

    python 06_crispr_library.py
"""
import pandas as pd

DATA = "data/essential_gene_gRNA_categorized.csv"


def main():
    g = pd.read_csv(DATA)
    print(f"essential genes targeted: {len(g)}")

    # guides per functional category
    print("\ngenes per category:")
    print(g["category"].value_counts().to_string())

    # count the guides that were actually designed (a few genes are incomplete)
    n_top = g["gRNA_top"].notna().sum()
    n_bot = g["gRNA_bottom"].notna().sum()
    both = (g["gRNA_top"].notna() & g["gRNA_bottom"].notna()).sum()
    print(f"\ngenes with both strands covered: {both} / {len(g)}")
    print(f"guides designed: {n_top} top-strand + {n_bot} bottom-strand = {n_top + n_bot}")

    # quick sanity check on guide length and PAM (GG = the last two bases of NGG)
    complete = g.dropna(subset=["g1_len", "g1_PAM"])
    print(f"median guide length: {complete['g1_len'].median():.0f} nt")
    print(f"guides ending in a GG PAM: {(complete['g1_PAM'] == 'GG').sum()} / {len(complete)}")


if __name__ == "__main__":
    main()
