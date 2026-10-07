#!/bin/bash
# ResFinder 4.7.2 (conda env "resfinder"), ResFinder database 2.6.0 (commit eecf0aa2). Acquired genes only.
# VERIFIED: the per-genome command is read from the software_executions entry of the ResFinder JSON of the study run:
#   __main__.py -ifa highquality_fastas/<acc>.fna -o results/resfinder/<acc> -s "Acinetobacter baumannii" -db_res <resfinder_db> -acq -l 0.6 -t 0.9
#   i.e. acquired resistance only (PointFinder and DisinFinder were not run, hence the empty pointfinder_blast / disinfinder_blast
#   folders), minimum coverage 60 % (-l 0.6), minimum identity 90 % (-t 0.9). The batch loop itself is a reconstruction.
# Set DB_RES to the local clone of the ResFinder database (the study used a database folder shared with another project).
set -euo pipefail
DB_RES=${DB_RES:-resfinder_db}
mkdir -p results/resfinder
while read -r ACC; do
  [[ -d results/resfinder/${ACC} && -s results/resfinder/${ACC}/ResFinder_results_tab.txt ]] && continue
  run_resfinder.py -ifa highquality_fastas/${ACC}.fna -o results/resfinder/${ACC} \
        -s "Acinetobacter baumannii" -db_res "$DB_RES" -acq -l 0.6 -t 0.9
done < highquality_accessions.txt
