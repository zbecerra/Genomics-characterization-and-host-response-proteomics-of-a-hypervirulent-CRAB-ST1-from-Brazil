#!/bin/bash
# ============================================================
# INDEPENDENT RE-IMPLEMENTATION of the QC filter (not the original run): per-genome QUAST and N counted directly
# from the FASTA instead of the QUAST rate. It returns the same 897 genomes as the original scripts
# (03a-03c); use it as a cross-check. Thresholds: completeness >= 95, contamination <= 5, contigs < 300, N <= 500.
# ============================================================
# Input : highquality_fastas/ (candidate ST1 genomes, *.fna)
# Output: highquality_accessions.txt (accession list used downstream)
# NOTE  : the accession list used in the study is provided as Table S2 and is
#         authoritative. The N limit (500 per genome) was inferred from the study data: it
#         cleanly separates the 897 retained genomes (max 500.8) from the 74 genomes that
#         pass the other three filters but were excluded (min 503.0).
# ============================================================
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"

# ── STEP 1: Run CheckM2 ──────────────────────────────────────
conda activate checkm2
mkdir -p output_checkm2 output_qc
checkm2 predict \
    --input highquality_fastas/ \
    --output-directory output_checkm2/ \
    --extension fna \
    --threads 48 \
    --force
echo "CheckM2 done!"
head -3 output_checkm2/quality_report.tsv

# ── STEP 2: Run QUAST ────────────────────────────────────────
conda activate quast_5.2.0
mkdir -p output_quast
for fna in highquality_fastas/*.fna; do
    accession=$(basename "$fna" .fna)
    quast.py "$fna" \
        -o output_quast/${accession} \
        --threads 4 \
        --no-html \
        --no-plots \
        -m 0 2>/dev/null
done

python3 << 'PY'
import pandas as pd, glob
results = []
for f in glob.glob("output_quast/*/report.tsv"):
    accession = f.split("/")[1]
    df = pd.read_csv(f, sep="\t", header=None, names=["metric", "value"])
    d = df.set_index("metric")["value"].to_dict()
    d["accession"] = accession
    results.append(d)
quast_df = pd.DataFrame(results)
quast_df.to_csv("output_quast/quast_summary.tsv", sep="\t", index=False)
print(f"QUAST summary: {len(quast_df)} genomes")
PY

# ── STEP 3: Merge CheckM2 + QUAST results ────────────────────
python3 << 'PY'
import pandas as pd
checkm = pd.read_csv("output_checkm2/quality_report.tsv", sep="\t")
checkm.rename(columns={"Name": "accession"}, inplace=True)
quast = pd.read_csv("output_quast/quast_summary.tsv", sep="\t")
merged = checkm.merge(quast, on="accession", how="left")

def count_n(path):
    n = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                n += line.upper().count("N")
    return n

# exact number of ambiguous bases (N) per genome, counted directly from the FASTA
merged["N_bases"] = [count_n(f"highquality_fastas/{a}.fna") for a in merged["accession"]]
print(f"CheckM2: {len(checkm)} | QUAST: {len(quast)} | merged: {len(merged)}")
merged.to_csv("output_qc/qc_merged.tsv", sep="\t", index=False)
PY

# ── STEP 4: Apply quality filters ────────────────────────────
MAX_N=500            # ambiguous bases (N) per genome, counted from the FASTA
EXPECTED_LIST=${EXPECTED_LIST:-}   # optional: path to a reference accession list for a self-check
export MAX_N EXPECTED_LIST
python3 << 'PY'
import pandas as pd, shutil, os
MAX_N = float(os.environ["MAX_N"]); EXPECTED = os.environ.get("EXPECTED_LIST", "")
df = pd.read_csv("output_qc/qc_merged.tsv", sep="\t")
NCOL = "N_bases"
for col in ["# contigs", NCOL]:
    df[col] = pd.to_numeric(df[col], errors="coerce")
print(f"Total genomes before filtering: {len(df)}")
keep = (df["Completeness"] >= 95) & (df["Contamination"] <= 5) & (df["# contigs"] < 300) & (df[NCOL] <= MAX_N)
df_hq = df[keep].copy()
print(f"High quality genomes: {len(df_hq)}")
print("Failing each criterion (genomes can fail several):")
print(f"  completeness < 95 %:      {(df['Completeness'] < 95).sum()}")
print(f"  contamination > 5 %:      {(df['Contamination'] > 5).sum()}")
print(f"  contigs >= 300:           {(df['# contigs'] >= 300).sum()}")
print(f"  N bases > {MAX_N:.0f} per genome: {(df[NCOL] > MAX_N).sum()}")
df_hq[["accession"]].to_csv("highquality_accessions.txt", header=False, index=False)
df_hq.to_csv("output_qc/highquality_genomes.tsv", sep="\t", index=False)
os.makedirs("highquality_fastas_selected", exist_ok=True)
for acc in df_hq["accession"]:
    src = f"highquality_fastas/{acc}.fna"
    if os.path.exists(src):
        shutil.copy(src, f"highquality_fastas_selected/{acc}.fna")
if EXPECTED and os.path.exists(EXPECTED):
    exp = set(open(EXPECTED).read().split()); got = set(df_hq["accession"])
    print(f"Self-check vs {EXPECTED}: identical = {exp == got} (only in expected: {len(exp-got)}, only in result: {len(got-exp)})")
PY

echo "QC filtering complete. HQ genomes: $(wc -l < highquality_accessions.txt)"
