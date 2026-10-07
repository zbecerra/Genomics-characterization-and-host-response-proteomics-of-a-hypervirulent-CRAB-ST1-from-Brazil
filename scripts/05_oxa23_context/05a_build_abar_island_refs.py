#!/usr/bin/env python3
"""
05a_build_abar_island_refs.py   (helper of 05_abar_genome_screen.sh; no seqkit needed)
Builds the island reference set used for the genome-wide AbaR screen: one sequence per AbaR column.
The AbaR3 island is cut from the AB0057 chromosome (NC_011586.2) at 235,738-326,347 (1-based, inclusive);
the other references are used as downloaded (first record of each file).

Usage: 05a_build_abar_island_refs.py databases/abar results/abar_island_screen/AbaR_island_refs.fna
"""
import os
import sys

db, out = sys.argv[1], sys.argv[2]
# (name used in the AbaR matrix, file in databases/abar, optional 1-based inclusive cut)
REFS = [
    ("AbaR3/Tn6019", "AbaR3_NC_011586.fna", (235738, 326347)),
    ("AbaR4",        "AbaR4_HQ700358.fna",  None),
    ("AbaR4a",       "AbaR4a_JN129845.fna", None),
    ("AbaR4b",       "AbaR4b_JN129846.fna", None),
    ("AbaR4c",       "AbaR4c_JN129847.fna", None),
    ("AbaR25",       "AbaR25_JX481978.fna", None),
]

def first_record(path):
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    break
                header = line[1:].strip()
            else:
                seq.append(line.strip())
    return header, "".join(seq)

with open(out, "w") as o:
    for name, fname, cut in REFS:
        header, seq = first_record(os.path.join(db, fname))
        if cut:
            assert "NC_011586" in header, f"unexpected AbaR3 source: {header}"
            seq = seq[cut[0] - 1:cut[1]]
        o.write(f">{name}\n")
        for i in range(0, len(seq), 70):
            o.write(seq[i:i + 70] + "\n")
        print(f"{name:<14}{len(seq):>9} bp   from {fname}  [{header[:60]}]")
