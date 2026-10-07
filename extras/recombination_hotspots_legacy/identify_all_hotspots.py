#!/usr/bin/env python3
"""
identify_all_hotspots.py
For every recombination hotspot window (from recomb_density_by_window_v2.tsv),
extract the overlapping gene(s) from Panaroo's core_alignment_header.embl,
and flag any that match known AMR or virulence gene name patterns.
"""
import re

DENSITY = "recomb_density_by_window_v2.tsv"
EMBL = "../pangenome/panaroo_results/core_alignment_header.embl"
TOP_N = 20
MIN_EVENTS = 15  # only report windows with at least this many events

# Known AMR / virulence gene name prefixes (from your own TableS3/TableS5 panels)
AMR_PATTERNS = [
    r'^bla', r'^aac', r'^aph', r'^ant', r'^tet', r'^sul', r'^dfr', r'^mph',
    r'^msr', r'^cml', r'^floR', r'^fos', r'^adeC', r'^adeF', r'^adeG',
    r'^adeH', r'^merA', r'^merC', r'^merR', r'^qac', r'^gyrA', r'^parC',
    r'^pmrB', r'^armA', r'^catA', r'^cxpE'
]
VIR_PATTERNS = [
    r'^pil', r'^fim', r'^csu', r'^bap$', r'^ompA', r'^tsaP', r'^pga',
    r'^bfm', r'^abaI', r'^abaR$', r'^gsp', r'^tss', r'^tse4', r'^lpx',
    r'^lps', r'^gal', r'^wec', r'^neuC', r'^pse', r'^tagX', r'^bas',
    r'^bau', r'^bar', r'^entE', r'^hemO', r'^plc', r'^adeF', r'^adeG',
    r'^adeH'
]

def flag_gene(name):
    base = name.split('~~~')[0].split('_')[0]  # handle fused-gene labels like lpxP~~~lpxL
    for pat in AMR_PATTERNS:
        if re.match(pat, name, re.IGNORECASE) or re.match(pat, base, re.IGNORECASE):
            return "AMR"
    for pat in VIR_PATTERNS:
        if re.match(pat, name, re.IGNORECASE) or re.match(pat, base, re.IGNORECASE):
            return "VIRULENCE"
    return ""

# 1. Load density windows, sort, take top N
windows = []
with open(DENSITY) as f:
    next(f)
    for line in f:
        start, end, count = line.strip().split('\t')
        windows.append((int(start), int(end), int(count)))
windows.sort(key=lambda x: -x[2])
top_windows = [w for w in windows[:TOP_N] if w[2] >= MIN_EVENTS]

# 2. Load all EMBL features into memory
features = []
with open(EMBL) as f:
    lines = f.readlines()
i = 0
while i < len(lines):
    if lines[i].startswith("FT   feature"):
        parts = lines[i].split()
        rng = parts[2]
        start, end = map(int, rng.split(".."))
        label_line = lines[i+1]
        label = label_line.split("/label=")[1].strip()
        features.append((start, end, label))
        i += 2
    else:
        i += 1

print(f"Loaded {len(features)} core-gene alignment features.\n")
print(f"{'Window':<25}{'Events':<8}{'Genes (flag)':<60}")
print("-" * 95)

for wstart, wend, count in top_windows:
    overlapping = [lbl for (fs, fe, lbl) in features if fe >= wstart and fs <= wend]
    flagged = []
    for g in overlapping:
        f = flag_gene(g)
        flagged.append(f"{g}[{f}]" if f else g)
    gene_str = ", ".join(flagged)
    print(f"{wstart:,}-{wend:,}   {count:<8}{gene_str}")
