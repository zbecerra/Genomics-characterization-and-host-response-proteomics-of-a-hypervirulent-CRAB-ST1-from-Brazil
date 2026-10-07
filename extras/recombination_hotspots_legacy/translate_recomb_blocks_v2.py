#!/usr/bin/env python3
"""
translate_recomb_blocks_v2.py
Translates ClonalFrameML recombination blocks (SNP-alignment column space)
into REAL reference-genome coordinates using the Gubbins
summary_of_snp_distribution.vcf POS column as the true position lookup
(record N's POS = genome coordinate of SNP-alignment column N).
"""
import csv
from collections import defaultdict

VCF = "../phylogeny/acb_ST1.summary_of_snp_distribution.vcf"
BLOCKS = "acb_ST1_clonalframe.importation_status.txt"
BIN_SIZE = 20000  # 20kb windows, appropriate for a multi-Mb genome
OUT_TRANSLATED = "recomb_blocks_genome_coords_v2.tsv"
OUT_DENSITY = "recomb_density_by_window_v2.tsv"

print("Loading real genome positions from VCF...")
lookup = {}
with open(VCF) as f:
    idx = 1  # 1-based, matching ClonalFrameML's Beg/End numbering
    for line in f:
        if line.startswith("#"):
            continue
        pos = int(line.split("\t")[1])
        lookup[idx] = pos
        idx += 1

print(f"Loaded {len(lookup)} real genome positions (range: {min(lookup.values())}-{max(lookup.values())})")

print("Translating blocks...")
translated = []
skipped = 0
with open(BLOCKS) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        node = row["Node"]
        try:
            beg_aln = int(row["Beg"])
            end_aln = int(row["End"])
        except ValueError:
            skipped += 1
            continue
        beg_gen = lookup.get(beg_aln)
        end_gen = lookup.get(end_aln)
        if beg_gen is None or end_gen is None:
            skipped += 1
            continue
        lo, hi = min(beg_gen, end_gen), max(beg_gen, end_gen)
        translated.append((node, lo, hi))

print(f"Translated {len(translated)} blocks ({skipped} skipped).")

with open(OUT_TRANSLATED, "w") as f:
    f.write("Node\tGenome_Start\tGenome_End\n")
    for node, lo, hi in translated:
        f.write(f"{node}\t{lo}\t{hi}\n")
print(f"Wrote {OUT_TRANSLATED}")

if translated:
    max_pos = max(hi for _, _, hi in translated)
    n_bins = (max_pos // BIN_SIZE) + 1
    density = defaultdict(int)
    for node, lo, hi in translated:
        start_bin = lo // BIN_SIZE
        end_bin = hi // BIN_SIZE
        for b in range(start_bin, end_bin + 1):
            density[b] += 1

    with open(OUT_DENSITY, "w") as f:
        f.write("Window_Start\tWindow_End\tEvent_Count\n")
        for b in range(n_bins):
            f.write(f"{b*BIN_SIZE}\t{(b+1)*BIN_SIZE}\t{density.get(b,0)}\n")
    print(f"Wrote {OUT_DENSITY} ({n_bins} windows, bin size {BIN_SIZE}bp)")

    top = sorted(density.items(), key=lambda x: -x[1])[:20]
    print("\nTop 20 recombination hotspot windows (genome-wide, real coordinates):")
    for b, count in top:
        print(f"  {b*BIN_SIZE:,}-{(b+1)*BIN_SIZE:,}: {count} events")
else:
    print("No blocks translated - check file paths/formats.")
