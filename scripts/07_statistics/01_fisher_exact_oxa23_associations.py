#!/usr/bin/env python3
"""
01_fisher_exact_oxa23_associations.py
Fisher's exact test (SciPy v1.18.1) of genome-level association between
blaOXA-23 carriage and mobile genetic elements / resistance-island subtypes
across all 897 ST1 genomes (2x2 tables; Fisher chosen over chi-square because
several cells are sparse or zero).

Inputs: Supplementary workbook with
  - 'TableS4. AMR matrix.'  (row 2 = headers incl. 'blaOXA-23'; col A = accession)
  - 'TableS5. MGE matrix.'  (row 2 = headers incl. ISAba1, AbaR3/Tn6019, AbaR4, AbaR4a/b/c, AbaR25)
  - optional 2nd argument: results/abar_island_screen/abar_presence.tsv (output of scripts/05_oxa23_context/05_abar_genome_screen.sh);
    when given, the AbaR calls are taken from it instead of the workbook, and ISAba1 still comes from the workbook.
Usage: python 01_fisher_exact_oxa23_associations.py Supplementary_tables.xlsx [abar_presence.tsv]
"""
import sys
import openpyxl
from scipy.stats import fisher_exact

xlsx = sys.argv[1]
wb = openpyxl.load_workbook(xlsx, data_only=True)

ws3 = wb['TableS4. AMR matrix.']
oxa_col = next(c for c in range(1, ws3.max_column + 1) if ws3.cell(row=2, column=c).value == 'blaOXA-23')
oxa23 = {}
for r in range(3, ws3.max_row + 1):
    acc = ws3.cell(row=r, column=1).value
    if acc:
        v = ws3.cell(row=r, column=oxa_col).value
        oxa23[acc] = 0 if v in (None, '') else int(float(v))

ws4 = wb['TableS5. MGE matrix.']
targets = ['ISAba1', 'AbaR3/Tn6019', 'AbaR4', 'AbaR4a', 'AbaR25']
cols = {ws4.cell(row=2, column=c).value: c for c in range(1, ws4.max_column + 1)
        if ws4.cell(row=2, column=c).value in targets}
mge = {}
for r in range(3, ws4.max_row + 1):
    acc = ws4.cell(row=r, column=1).value
    if acc:
        mge[acc] = {n: (0 if ws4.cell(row=r, column=c).value in (None, '')
                        else int(float(ws4.cell(row=r, column=c).value))) for n, c in cols.items()}

if len(sys.argv) > 2:
    with open(sys.argv[2]) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            f = line.rstrip('\n').split('\t')
            mge.setdefault(f[0], {})
            for n, v in zip(hdr[1:], f[1:]):
                mge[f[0]][n] = int(v)

# AbaR4, AbaR4a and AbaR25 overlap almost completely: also test them as one AbaR4-type group
if len(sys.argv) > 2:
    for acc in mge:
        mge[acc]['AbaR4-type'] = int(any(mge[acc].get(k, 0) for k in ('AbaR4', 'AbaR4a', 'AbaR25')))
    targets.append('AbaR4-type')

common = set(oxa23) & set(mge)
print(f"Genomes analysed: {len(common)}")
print(f"{'MGE':<15}{'O+/M+':>8}{'O+/M-':>8}{'O-/M+':>8}{'O-/M-':>8}{'OR':>10}{'p-value':>12}")
for name in targets:
    a = b = c = d = 0
    for acc in common:
        o, m = oxa23[acc], mge[acc][name]
        a += (o == 1 and m == 1); b += (o == 1 and m == 0)
        c += (o == 0 and m == 1); d += (o == 0 and m == 0)
    try:
        odds, p = fisher_exact([[a, b], [c, d]])
    except Exception:
        odds, p = float('nan'), float('nan')
    print(f"{name:<15}{a:>8}{b:>8}{c:>8}{d:>8}{odds:>10.2f}{p:>12.2e}")
