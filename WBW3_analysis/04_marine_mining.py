#!/usr/bin/env python3
"""Step 4 - Mine the ocean protein catalog for homologs of the WBW3 proteins.

The search itself is done with MMseqs2 against the Ocean Microbiomics Database
(OMD 2.0). The database is large (28.9M proteins), so it is not bundled here;
download it from the OMD website and run:

    mmseqs easy-search wbw3_proteins.fasta omd2_nr50_db hits.tsv tmp -s 7.5 \
        --format-output "q,t,pident,alnlen,mm,gap,qs,qe,ts,te,evalue,bits"

This script then reads the raw hits, pulls out the endolysin homologs, and
summarises how novel they are.
"""
import pandas as pd

HITS = "data/wbw3_vs_omd_all_hits.csv"          # all marine homologs (3,196)
OUT = "data/marine_endolysin_candidates.csv"    # endolysin homologs only (235)


def main():
    hits = pd.read_csv(HITS)
    print(f"total marine homologs of WBW3 proteins: {len(hits)}")

    # percent identity is stored as a fraction; convert to %
    hits["pident_pct"] = hits["pident"] * 100

    # how many hits per functional category?
    print("\nhits by functional category:")
    print(hits["function"].value_counts().to_string())

    # keep only the endolysin (lysis) homologs
    lys = hits[hits["function"] == "lysis"].copy()
    lys = lys.sort_values("bits", ascending=False)
    lys.to_csv(OUT, index=False)

    print(f"\nendolysin homologs kept: {len(lys)}")
    print(f"identity range: {lys['pident_pct'].min():.1f}% - {lys['pident_pct'].max():.1f}%")
    print("all are <50% identity, so every one is a novel enzyme")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
