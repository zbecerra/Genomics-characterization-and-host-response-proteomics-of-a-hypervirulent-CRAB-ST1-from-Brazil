#!/bin/bash
# IQ-TREE command: verbatim from results/clonalframeml_corrected/core_897.log (IQ-TREE 3.1.2).
# ClonalFrameML 1.20: settings read from cfml_897_newtree.log (emsim = 100, show_progress = true, 38 threads,
# 897 sequences x 2,136,645 sites, 58,420 variable sites analysed). The exact flag spelling of the thread option
# is not recorded (see `ClonalFrameML -h`); an earlier 32-thread attempt (job 75939) was cancelled on 2 Oct 2026.
# Must run on the FULL alignment (invariant sites retained), never on the SNP-only alignment.
# Study run (5 Oct 2026): R/theta 0.618892, 1/delta 0.00202071 (delta 494.9 bp), nu 0.0237018 -> r/m 7.26.
set -euo pipefail
mkdir -p results/clonalframeml_corrected && cd results/clonalframeml_corrected
ALN=../pangenome/panaroo_results/core_gene_alignment_filtered_clean.aln
iqtree3 -s $ALN -m GTR+G -T 36 -pre core_897 -fast
ClonalFrameML core_897.treefile $ALN ST1_cfml_897_newtree -emsim 100 -show_progress true
# add the thread option here (38 threads were used)
# r/m = R/theta x delta x nu (global values in ST1_cfml_897_newtree.em.txt)
