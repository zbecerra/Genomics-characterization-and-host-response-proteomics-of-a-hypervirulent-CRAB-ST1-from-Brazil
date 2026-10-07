#!/bin/bash
# Run CheckM2 and QUAST on the 1,256 ST1 candidate genomes (kept_copy_ST1/*.fna, after deduplication).
# QUAST 5.2.0 (conda env quast_5.2.0): command read from output_assacineto/quast.log - ONE call on all genomes,
#   38 threads, default options (the report still contains "# contigs (>= 0 bp)" and "# N's per 100 kbp").
# CheckM2 1.1.0: version, 1,256 bins and 20 threads read from checkm2.log (29 Jun 2026); the call itself is
#   RECONSTRUCTED (input directory/extension not recorded in the log).
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"

conda activate quast_5.2.0
quast.py kept_copy_ST1/*.fna --output-dir output_assacineto/ --threads 38
# -> output_assacineto/transposed_report.tsv (used by the merge step)

conda activate checkm2
checkm2 predict --input kept_copy_ST1 -x fna --threads 20 \
        --output-directory output_assacineto/checkm2_output
# -> output_assacineto/checkm2_output/quality_report.tsv
