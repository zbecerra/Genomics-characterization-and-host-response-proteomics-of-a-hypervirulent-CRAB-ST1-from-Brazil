#!/usr/bin/env bash
# identify_oxa23_context.sh
# Pulls representative sequences for the top clusters, downloads the AB0057
# reference (CP001182.2, has a confirmed intact Tn2006/AbaR4), and BLASTs
# to identify which structure each cluster corresponds to.
set -euo pipefail

OUTDIR="oxa23_context"
TOP_CLUSTERS="0 2 42 44"   # dominant cluster + top 3 secondary clusters

echo ">> Step 1: Extract representative header for each target cluster"
> "${OUTDIR}/top_cluster_reps.txt"
for cl in $TOP_CLUSTERS; do
    awk -v cl="$cl" '
        /^>Cluster/ { in_cluster = ($2 == cl) ? 1 : 0; next }
        in_cluster && /\*$/ {
            # representative line looks like: 0  1234nt, >HEADER... *
            match($0, />[^.]+\.\.\./)
            h = substr($0, RSTART+1, RLENGTH-4)
            print h
        }
    ' "${OUTDIR}/oxa23_contexts_clustered.fasta.clstr" >> "${OUTDIR}/top_cluster_reps.txt"
done
echo "Representative headers:"
cat "${OUTDIR}/top_cluster_reps.txt"

echo ""
echo ">> Step 2: Extract those sequences into one FASTA"
> "${OUTDIR}/top_cluster_reps.fasta"
while read -r header; do
    seqkit grep -n -p "$header" "${OUTDIR}/oxa23_contexts_clustered.fasta" >> "${OUTDIR}/top_cluster_reps.fasta"
done < "${OUTDIR}/top_cluster_reps.txt"
grep -c ">" "${OUTDIR}/top_cluster_reps.fasta"

echo ""
echo ">> Step 3: Download AB0057 reference chromosome (CP001182.2 - confirmed Tn2006/AbaR4)"
efetch -db nuccore -id CP001182.2 -format fasta > "${OUTDIR}/AB0057_reference.fasta"
grep ">" "${OUTDIR}/AB0057_reference.fasta"

echo ""
echo ">> Step 4: BLAST each cluster representative against the reference"
blastn -query "${OUTDIR}/top_cluster_reps.fasta" \
       -subject "${OUTDIR}/AB0057_reference.fasta" \
       -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
       -out "${OUTDIR}/top_clusters_vs_AB0057.tsv"

echo "Results (query=your cluster rep, sseqid=AB0057, pident=%identity, qcovs=%query coverage):"
column -t "${OUTDIR}/top_clusters_vs_AB0057.tsv"
