#!/bin/bash
# Insertion-sequence search used for the named-IS columns of the MGE matrix (Table S5), plus the PlasmidFinder run (not used in S5).
# COMMAND TEMPLATE: the loop only logged "Done: <acc>" (results/mge/is_tn_blast, results/plasmids/plasmidfinder_run.log);
# the exact options were not recorded. What the outputs and the table show (verified):
#  - Insertion sequences: BLASTN of the genome FASTA against a nucleotide database built from ISFINDER ONLY (IS.fna, 5,970 sequences,
#    downloaded 2026-07-28). The file called IS_Tn_combined.fna is byte-identical to IS.fna (TnCentral was NOT added).
#    Columns: qseqid sseqid pident length qstart qend sstart send evalue bitscore, plus accession, family (IS name).
#    BLAST+ 2.5.0+ (conda env "mge"). A genome is scored 1 for an IS name in Table S5 when it has a hit of >= 80 % identity
#    (this rule reproduces the table for all 897 genomes).
#  - PlasmidFinder 2.1.6 (BLAST) was run with the Enterobacteriaceae and Gram-positive databases only: hits in 8/897 genomes (non-Acinetobacter
#    replicons). It is NOT used for Table S5; replicon/relaxase/MPF information in the table comes from MOB-suite.
set -euo pipefail
mkdir -p results/mge/is_tn_blast results/plasmids/plasmidfinder_results
while read -r ACC; do
  blastn -query highquality_fastas/${ACC}.fna -db databases/mge_combined/IS_Tn_db \
          -out results/mge/is_tn_blast/${ACC}.tsv \
          -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"    # e-value/identity options: not recorded
  plasmidfinder.py -i highquality_fastas/${ACC}.fna -o results/plasmids/plasmidfinder_results/${ACC} -p <plasmidfinder_db> -x  # options: assumed
done < highquality_accessions.txt
