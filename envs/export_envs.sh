#!/bin/bash
# Run ON THE CLUSTER where the conda environments live; writes one portable file per environment into envs/.
# --no-builds keeps the files installable on other systems. Environments that do not exist are skipped.
set -uo pipefail
mkdir -p envs
for e in mlst download panaroo snpsites quast_5.2.0 checkm2 gubbins clonalframe iqtree amrfinder rgi card_new resfinder abricate kaptive mobsuite mge integron islandpath plasmidfinder alignment phylo; do
  if conda env list | awk '{print $1}' | grep -qx "$e"; then conda env export -n "$e" --no-builds > "envs/$e.yml" && echo "wrote envs/$e.yml"; else echo "skip: no environment $e"; fi
done
# Bakta lives in the system (spack) conda installation; export it by path (adjust the path to your installation):
# conda env export -p /home/public/spack/opt/spack/linux-rocky9-zen3/gcc-12.2.0/miniconda3-24.3.0-afr7dkwcon7ppaf4mcdgnje2inigdmrt/envs/bakta --no-builds > envs/bakta.yml
