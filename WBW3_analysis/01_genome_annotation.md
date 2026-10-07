# Step 1 — Genome re-annotation

The original WBW3 annotation (RAST/BLAST, 76 ORFs) was redone with phage-specific
tools to catch genes that plain sequence similarity misses.

## 1.1 Re-orient the genome

The circular genome is re-oriented to start at the terminase large subunit
(the conventional starting point for phage genomes).

```bash
dnaapler --input WBW3.fasta --output WBW3_reoriented.fasta --mode phage
```

## 1.2 Pharokka (gene calling + first-pass annotation)

Pharokka runs PHANOTATE for gene calling and screens against PHROGs, CARD and
VFDB (for resistance/virulence genes), plus tRNAscan-SE and MinCED.

```bash
pharokka --genome WBW3_reoriented.fasta --outdir pharokka_out --prefix WBW3
```

This gave 83 protein-coding genes (up from the original 76).

## 1.3 Phold (structure-informed annotation)

Phold uses a protein language model (ProstT5) to predict 3Di structure and then
Foldseek to search a structure database. This is what resolved most of the
"hypothetical" proteins that sequence search alone could not annotate.

```bash
phold --input pharokka_out/WBW3.gbk --outdir phold_out --prefix WBW3
```

## 1.4 Merge and tidy the annotation

The Pharokka + Phold calls were merged into one table. The final per-gene table
is `data/WBW3_final_annotation_pharokka_phold.csv` (83 CDS, with functional
category, product, and a confidence label for each call).

Result: 45 of 83 CDS (54%) got a functional assignment, versus ~40% with the
original RAST/BLAST annotation. Eight previously unknown ORFs were resolved by
the structure-based search (portal protein, DNA helicase, nuclease, nucleotide
kinase, glycosyltransferase, etc.). No resistance or virulence genes were found.
