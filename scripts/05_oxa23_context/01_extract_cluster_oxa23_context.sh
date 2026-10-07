#!/usr/bin/env bash
# oxa23_context_extraction.sh (v2 - uses Bakta coordinates)
#
# AMRFinder was run in protein mode (no contig/start/stop columns), so this
# version joins AMRFinder's "Protein id" against Bakta's "Locus Tag" to
# recover genomic coordinates, then extracts and clusters flanking context
# around every blaOXA-23 hit.
set -euo pipefail

AMRFINDER_DIR="results/amr/amrfinder_results"
BAKTA_DIR="bakta_output"
FLANK=5000
OUTDIR="oxa23_context"

mkdir -p "${OUTDIR}"

echo ">> Step 1: Find blaOXA-23 hits (Protein id) per genome"
> "${OUTDIR}/oxa23_hits_proteinid.tsv"
for f in "${AMRFINDER_DIR}"/*.tsv; do
    acc=$(basename "$f" .tsv)
    awk -F'\t' -v acc="$acc" '
        tolower($2) ~ /^blaoxa$/ && tolower($3) ~ /oxa-23/ {print acc"\t"$1}
        tolower($2) ~ /oxa-23/ {print acc"\t"$1}
    ' "$f" >> "${OUTDIR}/oxa23_hits_proteinid.tsv"
done
sort -u -o "${OUTDIR}/oxa23_hits_proteinid.tsv" "${OUTDIR}/oxa23_hits_proteinid.tsv"
n_hits=$(wc -l < "${OUTDIR}/oxa23_hits_proteinid.tsv")
echo "Found ${n_hits} blaOXA-23 hits."
head -5 "${OUTDIR}/oxa23_hits_proteinid.tsv"

echo ""
echo ">> Step 2: Join against Bakta TSV to get contig/start/stop"
> "${OUTDIR}/oxa23_hits_coords.tsv"
while IFS=$'\t' read -r acc protid; do
    bakta_tsv="${BAKTA_DIR}/${acc}/${acc}.tsv"
    if [[ ! -f "$bakta_tsv" ]]; then
        echo "  !! Missing Bakta TSV for ${acc}, skipping"
        continue
    fi
    # Bakta columns (after stripping # comments): Sequence_Id Type Start Stop Strand Locus_Tag Gene Product DbXrefs
    grep -v "^#" "$bakta_tsv" | awk -F'\t' -v acc="$acc" -v pid="$protid" \
        '$6==pid {print acc"\t"$1"\t"$3"\t"$4}' >> "${OUTDIR}/oxa23_hits_coords.tsv"
done < "${OUTDIR}/oxa23_hits_proteinid.tsv"

n_coords=$(wc -l < "${OUTDIR}/oxa23_hits_coords.tsv")
echo "Resolved coordinates for ${n_coords}/${n_hits} hits."
head -5 "${OUTDIR}/oxa23_hits_coords.tsv"

echo ""
echo ">> Step 3: Extract ${FLANK} bp flanking each hit (from Bakta's own .fna, contig naming matches)"
> "${OUTDIR}/oxa23_contexts.fasta"
while IFS=$'\t' read -r acc contig start stop; do
    fasta="${BAKTA_DIR}/${acc}/${acc}.fna"
    if [[ ! -f "$fasta" ]]; then
        echo "  !! Missing Bakta fna for ${acc}, skipping"
        continue
    fi
    win_start=$(( start - FLANK > 0 ? start - FLANK : 1 ))
    win_end=$(( stop + FLANK ))
    seqkit faidx "$fasta" "${contig}:${win_start}-${win_end}" 2>/dev/null \
        | seqkit replace -p "^(.+)$" -r "${acc}_${contig}_${win_start}-${win_end}" \
        >> "${OUTDIR}/oxa23_contexts.fasta"
done < "${OUTDIR}/oxa23_hits_coords.tsv"

n_extracted=$(grep -c ">" "${OUTDIR}/oxa23_contexts.fasta" || echo 0)
echo "Extracted ${n_extracted} flanking-region sequences."

echo ""
echo ">> Step 4: Cluster contexts by similarity (95% identity, 90% coverage)"
cd-hit-est -i "${OUTDIR}/oxa23_contexts.fasta" \
           -o "${OUTDIR}/oxa23_contexts_clustered.fasta" \
           -c 0.95 -aS 0.90 -d 0 -T 4 -M 4000

n_clusters=$(grep -c ">" "${OUTDIR}/oxa23_contexts_clustered.fasta")
echo "Collapsed into ${n_clusters} representative genetic contexts."

echo ""
echo ">> Step 5: Cluster size summary"
awk '/^>Cluster/{if(name!="") print name"\t"count; name=$0; count=0; next} {count++} END{print name"\t"count}' \
    "${OUTDIR}/oxa23_contexts_clustered.fasta.clstr" | sort -t$'\t' -k2 -rn > "${OUTDIR}/cluster_sizes.tsv"
cat "${OUTDIR}/cluster_sizes.tsv"

echo ""
echo ">> Next (optional, run separately):"
echo "   mafft ${OUTDIR}/oxa23_contexts_clustered.fasta > ${OUTDIR}/oxa23_contexts_aligned.fasta"
echo "   blastn -query ${OUTDIR}/oxa23_contexts_clustered.fasta -subject tn_references.fasta -outfmt 6"
