#!/bin/bash
# RECONSTRUCTED (original command not recorded). Panaroo v1.7.0, MAFFT v7.526 (conda env "panaroo").
# Verified from the outputs: 897 genomes, 10,727 genes, 2,246 genes in >=99 % of strains, and a core alignment
# built with core_threshold 0.98 (2,519 genes) per alignment_resume_state.json.
# Clean mode and thread count of the original run are NOT recorded.
# Run from the project root. Input: Bakta GFF3 of the 897 genomes in highquality_accessions.txt.
set -euo pipefail
mkdir -p results/pangenome
GFFS=$(while read -r a; do echo "bakta_output/${a}/${a}.gff3"; done < highquality_accessions.txt)
panaroo -i $GFFS -o results/pangenome/panaroo_results \
        --clean-mode strict -a core --aligner mafft --core_threshold 0.98 -t 16
# Panaroo writes core_gene_alignment.aln and the filtered core_gene_alignment_filtered.aln (used downstream).
# The _clean file differs from the filtered one by 5 ambiguity characters (w, s, m, m, k -> N):
cd results/pangenome/panaroo_results
python3 - <<'PY'
import re
src, dst = "core_gene_alignment_filtered.aln", "core_gene_alignment_filtered_clean.aln"
with open(src) as f, open(dst, "w") as o:
    for line in f:
        o.write(line if line.startswith(">") else re.sub(r"[^ACGTNacgtn\-]", "N", line))
PY
