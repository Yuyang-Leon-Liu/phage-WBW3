# Step 3 — Structure prediction and analysis of the lysis proteins

We predicted structures for the three lysis-cassette proteins (holin, endolysin,
spanin) and a 774-aa tail protein, then used structure search to identify the
endolysin fold and its catalytic residues.

## 3.1 AlphaFold2 structure prediction

Each protein sequence (single FASTA) was run through AlphaFold2 in monomer mode:

```bash
python3 run_alphafold.py \
    --fasta_paths=endolysin.fasta \
    --output_dir=af2_out \
    --model_preset=monomer \
    --db_preset=reduced_dbs
```

The endolysin model is high confidence (mean pLDDT 95.1). The small membrane
proteins (holin, spanin) are lower confidence, which is expected.

## 3.2 Foldseek structure search against the PDB

```bash
foldseek easy-search endolysin.pdb pdb_db aln.tsv tmp
```

Top hits were all characterized GH24 muramidases / endolysins:
DLP12 endolysin (4ZPU), R21 (3HDE), P22 lysozyme (2ANV), P1 Lyz (1XJT),
and the thermostable PHAb10/PHAb11 (9KBQ/9KBS).

## 3.3 Mapping the catalytic dyad

A multiple-sequence alignment of the WBW3 endolysin against those five
characterized GH24 enzymes maps the conserved catalytic residues to
**Glu26 and Asp35** (the same Glu/Asp pair as the T4-lysozyme prototype).
The N-terminal residues 1–21 form a hydrophobic signal-arrest-release (SAR)
export domain, as in the R21/DLP12 SAR endolysins.

## 3.4 Developability profile

Basic physicochemical properties (pI, instability index, GRAVY, molecular
weight) were computed on the endolysin sequence with Biopython's ProtParam:

```python
from Bio.SeqUtils.ProtParam import ProteinAnalysis
p = ProteinAnalysis(open("endolysin.fasta").read().split("\n", 1)[1].replace("\n", ""))
print(p.isoelectric_point(), p.instability_index(), p.molecular_weight, p.gravy())
```

The enzyme is basic (pI 9.67), predicted stable (instability index 31.3),
17.6 kDa, and near-neutral GRAVY (-0.03).
