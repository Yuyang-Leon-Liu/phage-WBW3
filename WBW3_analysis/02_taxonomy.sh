#!/usr/bin/env bash
# Step 2 — Taxonomic placement of WBW3 within Drexlerviridae.
# Three independent lines of evidence: whole-genome ANI, a terminase phylogeny,
# and a gene-sharing network. Reference genomes are from INPHARED.
set -euo pipefail

# --- inputs ---------------------------------------------------------------
WBW3=WBW3_reoriented.fasta          # re-oriented genome from step 1
INPHARED=inphared_genomes.fasta     # 34,076 phage genomes (14 Apr 2025 release)
                                    # download: https://github.com/RyanCook94/inphared

# --- 2.1 whole-genome ANI (skani) ----------------------------------------
# A new species is <95% ANI. skani is fast enough for an all-vs-all scan.
skani dist "$WBW3" "$INPHARED" -o wbw3_vs_inphared_ani.tsv
# -> max ANI was 88.8% (to Veterinaerplatzvirus vB_EcoS-2004IV), i.e. a new species.

# --- 2.2 gene-sharing network (MMseqs2 clustering) ------------------------
# Cluster all proteins at 30% identity / 50% coverage, then count how many
# clusters each reference genome shares with WBW3 (a vConTACT2-style approach).
mmseqs createdb "$WBW3" "$INPHARED" all_proteins_db
mmseqs linclust all_proteins_db clusters_db tmp \
    --min-seq-id 0.3 -c 0.5
mmseqs createtsv all_proteins_db clusters_db protein_clusters.tsv
# -> 40,403 proteins fell into 2,013 clusters; WBW3 shares the most (50)
#    with Escherichia phage vB_EcoS_CEB_EC3a.

# --- 2.3 terminase large-subunit phylogeny --------------------------------
# Align terminase proteins with MAFFT, then build a maximum-likelihood tree.
mafft --auto terminase_proteins.fasta > terminase_aligned.fasta
iqtree2 -s terminase_aligned.fasta -m MFP -B 1000 -T AUTO --prefix terl_tree
# -> WBW3 sits inside the Veterinaerplatzvirus clade (Braunvirinae).

echo "Done. Check wbw3_vs_inphared_ani.tsv, protein_clusters.tsv and terl_tree.treefile"
