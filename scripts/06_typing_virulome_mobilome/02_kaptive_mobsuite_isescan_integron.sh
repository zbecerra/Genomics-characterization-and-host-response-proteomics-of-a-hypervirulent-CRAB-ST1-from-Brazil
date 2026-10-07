#!/bin/bash
# COMMAND TEMPLATES - versions verified (Kaptive 2.0.6, MOB-suite 3.1.9, ISEScan 1.7.3 in conda env "mge",
# IntegronFinder 2.0.6, IslandPath-DIMOB 1.0.6 in env "islandpath"); the original commands were not recorded.
# Evidence from the logs: ISEScan ran per genome on highquality_fastas/<acc>.fna into results/isescan/<acc> with
# FragGeneScan -p 2 (so --nthread 2); IslandPath-DIMOB (Dimob.pl) ran on the Bakta .gbff files (bakta_output/<acc>/<acc>.gbff). Fill in the Kaptive locus database paths before use.
# IslandPath-DIMOB v1.0.6: https://github.com/brinkmanlab/islandpath
set -euo pipefail
KLOCUS="/path/to/Acinetobacter_baumannii_k_locus_primary_reference.gbk"
OCLOCUS="/path/to/Acinetobacter_baumannii_OC_locus_primary_reference.gbk"
mkdir -p results/kaptive results/mobsuite results/isescan results/integrons
# The study output has a per-genome FASTA, one JSON and one table per locus database, and no images, so --json and
# image suppression were probably used (check `kaptive.py --help`); the output prefixes below match the study files.
kaptive.py -a highquality_fastas/*.fna -k $KLOCUS  -o results/kaptive/acb_ST1_K_locus
kaptive.py -a highquality_fastas/*.fna -k $OCLOCUS -o results/kaptive/acb_ST1_OC_locus
while read -r ACC; do
  mob_recon -i highquality_fastas/${ACC}.fna -o results/mobsuite/${ACC} -n 8
  isescan.py --seqfile highquality_fastas/${ACC}.fna --output results/isescan/${ACC} --nthread 2
  integron_finder highquality_fastas/${ACC}.fna --outdir results/integrons --cpu 8   # creates results/integrons/Results_Integron_Finder_${ACC} (as in the study)
done < highquality_accessions.txt
# IslandPath-DIMOB (env islandpath): input = Bakta GenBank file; output file name is a placeholder
mkdir -p results/islandpath
while read -r ACC; do
  Dimob.pl bakta_output/${ACC}/${ACC}.gbff results/islandpath/${ACC}.dimob.txt
done < highquality_accessions.txt
