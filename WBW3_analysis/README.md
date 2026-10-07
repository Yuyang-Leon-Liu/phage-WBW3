# WBW3 analysis

Code and data for the phage genomics study of **vB_EcoD-WBW3**, a T1-like
*Escherichia coli* phage, and the mining of its endolysin against the ocean
protein catalog.

The pipeline has four parts: re-annotate the genome, place it taxonomically,
model the lysis proteins, and mine marine virome data for related enzymes.
Each step is in its own file, numbered in order.

## What's in here

| File | What it does |
|------|--------------|
| `01_genome_annotation.md` | Re-orient the genome and re-annotate it (Pharokka + Phold) |
| `02_taxonomy.sh` | ANI, gene-sharing network, and terminase phylogeny |
| `03_structure_analysis.md` | AlphaFold2 + Foldseek analysis of the lysis proteins |
| `04_marine_mining.py` | Pull the marine endolysin homologs out of the OMD search results |
| `05_make_figures.py` | Redraw the data-driven figures from the tables |
| `06_crispr_library.py` | Summarise the essential-gene CRISPR guide library |
| `data/` | The result tables (CSVs) used by the scripts |

## Running it

The Python scripts only need the tables in `data/`:

```bash
pip install -r requirements.txt
python 04_marine_mining.py
python 05_make_figures.py
python 06_crispr_library.py
```

Steps 1–3 use command-line bioinformatics tools (Pharokka, Phold, skani,
MMseqs2, MAFFT, IQ-TREE2, AlphaFold2, Foldseek). Those are run at the command
line as shown in the files; the large reference databases they need (INPHARED,
OMD 2.0, the Phold database) are downloaded from their own websites and are not
included here because of their size.

## Data sources

- **INPHARED** phage genomes — https://github.com/RyanCook94/inphared
- **Ocean Microbiomics Database (OMD 2.0)** — https://www.ocean-microbiome.com
- **Phold / Pharokka** — https://github.com/gbouras13/phold

## Notes

The WBW3 genome sequence is deposited in GenBank (see the paper's Data
Availability Statement for the accession number). If you use this code or data,
please cite the paper (see `CITATION.cff`).
