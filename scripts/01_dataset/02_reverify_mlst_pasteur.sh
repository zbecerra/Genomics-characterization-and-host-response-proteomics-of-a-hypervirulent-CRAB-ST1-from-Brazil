#!/usr/bin/env bash
# reverify_st1.sh
# Re-runs MLST typing on all genomes already in kept_ST1/ using the patched
# mlst 2.17.6 (BLAST+ version-check bug fixed), with full logging this time.
set -uo pipefail

FASTA_DIR="ACB_ST1_pipeline/kept_ST1"
OUTDIR="ACB_ST1_pipeline/logs/reverify"
mkdir -p "${OUTDIR}"

OUTFILE="${OUTDIR}/mlst_reverify_full.tsv"
echo -e "accession\tscheme\tST\talleles" > "${OUTFILE}"

n=0
total=$(ls "${FASTA_DIR}"/*.fna | wc -l)

for f in "${FASTA_DIR}"/*.fna; do
    n=$((n+1))
    acc=$(basename "$f" .fna)
    result=$(mlst --scheme abaumannii_2 --quiet "$f" 2>/dev/null)
    if [[ -z "$result" ]]; then
        echo -e "${acc}\tabaumannii_2\tFAILED\t-" >> "${OUTFILE}"
    else
        scheme=$(echo "$result" | awk -F'\t' '{print $2}')
        st=$(echo "$result" | awk -F'\t' '{print $3}')
        alleles=$(echo "$result" | cut -f4-)
        echo -e "${acc}\t${scheme}\t${st}\t${alleles}" >> "${OUTFILE}"
    fi
    if (( n % 100 == 0 )); then
        echo "Progress: ${n}/${total}" >&2
    fi
done

echo ""
echo "=== Done. Summary of ST calls across all ${total} genomes: ==="
awk -F'\t' 'NR>1{print $3}' "${OUTFILE}" | sort | uniq -c | sort -rn

echo ""
echo "=== Non-ST1 or failed genomes (should be investigated/excluded): ==="
awk -F'\t' 'NR>1 && $3!="1" {print}' "${OUTFILE}"
