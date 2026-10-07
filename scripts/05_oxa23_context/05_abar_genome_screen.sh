#!/usr/bin/env bash
# 05_abar_genome_screen.sh        STATUS: NEW ANALYSIS (replaces the earlier genome-wide AbaR screen; not yet run when written)
#
# Genome-wide presence/absence of the AbaR islands in the 897 genomes (AbaR3/Tn6019, AbaR4, AbaR4a, AbaR25; AbaR4b/4c are computed but not reported),
# used for the association tests (Step 9) and for the AbaR columns of Table S5.
#
# Why this replaces the earlier screen: that search used a database in which "AbaR3" was the WHOLE AB0057 chromosome
# (NC_011586.2, 4 Mb), so every A. baumannii genome matched it (897/897). Here every reference is the island itself and
# presence is a stated rule: >= 90 % nucleotide identity over >= 80 % of the reference length (union of aligned intervals).
#
# AbaR3 reference: the resistance-gene island of AB0057 (NC_011586.2 = Bakta contig_1 of GCF_000021245.2), positions
# 235,738-326,347 (Bakta annotation: tns/tni module, sul1, mer operon, tet(A), catA1, aph(3')-Ia, aadA1, aac(3)-Ia; ~90.6 kb).
# The other references are the characterised islands in databases/abar (AbaR4_HQ700358, AbaR4a_JN129845, AbaR4b_JN129846, AbaR4c_JN129847, AbaR25_JX481978).
# Note: AbaR4, AbaR4a, AbaR4b and AbaR4c share most of their sequence, so their calls are NOT mutually exclusive.
set -euo pipefail

DB_DIR="databases/abar"
OUT="results/abar_island_screen"
GENOMES="highquality_fastas"          # the 897 genome FASTA files (BLAST -subject mode: no per-genome database needed)
THREADS=8
MIN_ID=90
MIN_COV=0.80
mkdir -p "${OUT}"

echo ">> Step 1: Build the island reference set (one sequence per column of the AbaR matrix)"
# AbaR3 island cut from the AB0057 chromosome (235,738-326,347); the other references are used as downloaded.
python3 scripts/05_oxa23_context/05a_build_abar_island_refs.py "${DB_DIR}" "${OUT}/AbaR_island_refs.fna"

echo ">> Step 2: BLASTN of the references against every genome (reference = query, genome = subject)"
run_one() {
  acc="$1"
  blastn -query "${OUT}/AbaR_island_refs.fna" -subject "${GENOMES}/${acc}.fna" -perc_identity "${MIN_ID}" -evalue 1e-10 \
         -outfmt "6 qseqid sseqid pident length qstart qend sstart send qlen evalue bitscore" \
         | awk -v a="${acc}" 'BEGIN{OFS="\t"}{print $0, a}' > "${OUT}/hits_${acc}.tsv"
}
export -f run_one; export OUT GENOMES MIN_ID
xargs -P "${THREADS}" -I{} bash -c 'run_one {}' < highquality_accessions.txt     # the 897 accessions (not every FASTA in the folder)
cat "${OUT}"/hits_*.tsv > "${OUT}/abar_island_hits_all.tsv"
rm -f "${OUT}"/hits_*.tsv

echo ">> Step 3: Presence/absence (union of aligned intervals / reference length >= ${MIN_COV})"
python3 scripts/05_oxa23_context/05b_abar_presence.py \
  "${OUT}/abar_island_hits_all.tsv" highquality_accessions.txt "${MIN_COV}" \
  "${OUT}/abar_presence.tsv" "${OUT}/abar_fraction_covered.tsv"
