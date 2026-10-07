#!/usr/bin/env python3
"""
05b_abar_presence.py   (helper of 05_abar_genome_screen.sh)
Turns the BLASTN hits of the AbaR island references into a presence/absence table.

A reference is "present" in a genome when the union of its aligned intervals (query coordinates,
overlapping hits merged; hits already filtered at >= 90 % identity) covers >= MIN_COV of its length.

Usage: 05b_abar_presence.py hits.tsv accessions.txt MIN_COV presence.tsv fraction.tsv
hits.tsv columns: qseqid sseqid pident length qstart qend sstart send qlen evalue bitscore accession
"""
import sys
from collections import defaultdict

hits_f, acc_f, min_cov, out_presence, out_frac = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4], sys.argv[5]
REFS = ["AbaR3/Tn6019", "AbaR4", "AbaR4a", "AbaR4b", "AbaR4c", "AbaR25"]

accs = [l.strip() for l in open(acc_f) if l.strip()]
iv = defaultdict(list)           # (accession, ref) -> [(start, end)]
qlen = {}
for line in open(hits_f):
    f = line.rstrip("\n").split("\t")
    ref, qs, qe, ql, acc = f[0], int(f[4]), int(f[5]), int(f[8]), f[-1]
    qlen[ref] = ql
    iv[(acc, ref)].append((min(qs, qe), max(qs, qe)))

def covered(intervals):
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total

with open(out_presence, "w") as p, open(out_frac, "w") as q:
    p.write("accession\t" + "\t".join(REFS) + "\n")
    q.write("accession\t" + "\t".join(REFS) + "\n")
    for acc in accs:
        fr = [covered(iv.get((acc, r), [])) / qlen[r] if r in qlen else 0.0 for r in REFS]
        p.write(acc + "\t" + "\t".join("1" if x >= min_cov else "0" for x in fr) + "\n")
        q.write(acc + "\t" + "\t".join(f"{x:.3f}" for x in fr) + "\n")
n = len(accs)
print("genomes:", n)
for j, r in enumerate(REFS):
    c = sum(1 for l in open(out_presence).read().splitlines()[1:] if l.split("\t")[j + 1] == "1")
    print(f"{r:<14} present in {c}/{n}")
