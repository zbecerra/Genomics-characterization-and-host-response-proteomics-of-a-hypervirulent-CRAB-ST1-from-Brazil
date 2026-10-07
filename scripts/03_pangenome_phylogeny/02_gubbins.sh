#!/bin/bash
# RECONSTRUCTED from results/gubbins/gubbins_clean_run.log and results/phylogeny/acb_ST1.log.
# Gubbins 3.4.3, RAxML 8.2.12 (GTRGAMMA), 20 threads, input = core_gene_alignment_filtered_clean.aln.
# Verified results of the study run: 83,487 of 100,744 SNPs (82.87 %) in recombinant blocks; pooled r/m 4.84;
# mean per-branch r/m 2.93 over 1,793 branches; 4,347 blocks on 353 branches.
set -euo pipefail
mkdir -p results/gubbins && cd results/gubbins
run_gubbins.py --prefix acb_ST1 --tree-builder raxml --threads 20 \
    ../pangenome/panaroo_results/core_gene_alignment_filtered_clean.aln
# mean per-branch r/m
awk -F',' 'NR>1{s+=$2;n++} END{print "mean r/m:",s/n,"(n="n" branches)"}' acb_ST1.per_branch_statistics.csv
