#!/usr/bin/env bash
# 01_build_acinetobacter_vfdb_and_screen.sh
# Virulence screening with ABRicate against an Acinetobacter-specific subset of
# VFDB Set B (nucleotide, VFDB_setB_nt.fas; http://www.mgc.ac.cn/VFs/).
#
# History / fix: the first custom databases (acb_vfdb_full, acb_vf) had every
# sequence truncated to exactly 60 bp, which silently changed results. The
# database was rebuilt from the raw VFDB file (966 full-length sequences).
# NOTE: "abricate --db vfdb" (ABRicate's bundled VFDB) was NOT used for the
# manuscript results; thresholds are ABRicate defaults (>=80 % identity, >=80 % coverage).
set -euo pipefail

VFDB_RAW="VFDB_setB_nt.fas"
ABRICATE_DB_DIR="${CONDA_PREFIX}/db"          # abricate --help shows the datadir
DB_NAME="acb_vfdb_rebuilt"
FASTA_DIR="ACB_ST1_pipeline/kept_ST1"          # use the 897 high-quality genomes
OUT_DIR="vfdb_rebuilt_results"

# 1) keep only entries whose header names Acinetobacter as source organism
awk '/^>/{p = ($0 ~ /Acinetobacter/)} p' "${VFDB_RAW}" > acb_vfdb_rebuilt.fasta
grep -c ">" acb_vfdb_rebuilt.fasta             # expect 966

# 2) register as an ABRicate database
mkdir -p "${ABRICATE_DB_DIR}/${DB_NAME}"
cp acb_vfdb_rebuilt.fasta "${ABRICATE_DB_DIR}/${DB_NAME}/sequences"
abricate --setupdb
abricate --list | grep "${DB_NAME}"            # should report 966 sequences

# 3) screen genomes (default thresholds)  -- loop reconstructed from output file naming
mkdir -p "${OUT_DIR}"
for f in "${FASTA_DIR}"/*.fna; do
    acc=$(basename "${f}" .fna)
    abricate --db "${DB_NAME}" "${f}" > "${OUT_DIR}/${acc}_${DB_NAME}.tsv"
done
