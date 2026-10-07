#!/bin/bash
# Recombination-filtered cgSNP tree and SNP-distance matrix (outputs in results/phylogeny/).
# VERIFIED (iqtree_run.logn, IQ-TREE 3.1.2): the IQ-TREE command below, run on the Gubbins recombination-filtered
#   polymorphic sites with the constant-site counts (A,C,G,T) in -fconst and 1000 UFBoot2 replicates.
# VERIFIED from the output files: acb_ST1_snpsites.{fasta,phylip,vcf} (5 Jul 2026) = 897 sequences x 15,282 SNP sites;
#   acb_ST1_snp_distance_matrix.tsv (17 Jul 2026) = 897 x 897, first cell "snp-dists 1.2.0"; its input is acb_ST1_snpsites.fasta
#   (VERIFIED: three pairwise distances recounted from that FASTA, 119 / 93 / 68, equal the matrix values).
# NOT RECORDED: the SNP-sites flags and its input file (assumed: the Gubbins filtered polymorphic sites) and how the -fconst counts were derived (Gubbins documents `snp-sites -C`).
set -euo pipefail
cd results/phylogeny      # holds acb_ST1.filtered_polymorphic_sites.fasta (Gubbins output) in the study layout
mkdir -p iqtree_results

# SNP-sites (conda env "snpsites") - RECONSTRUCTED flags
snp-sites -m -v -p -o acb_ST1_snpsites acb_ST1.filtered_polymorphic_sites.fasta
# SNP-dists 1.2.0 (same env) - input verified
snp-dists acb_ST1_snpsites.fasta > acb_ST1_snp_distance_matrix.tsv

# IQ-TREE 3.1.2 - verbatim from iqtree_run.logn
iqtree -s acb_ST1.filtered_polymorphic_sites.fasta -m GTR+G -fconst 619984,389653,445493,622792 \
       -bb 1000 -nt 38 -pre iqtree_results/acb_ST1
