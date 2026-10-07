#!/bin/bash
# AMRFinderPlus v4.2.7 (DB 2026-05-15.1): command taken from results/amr/amrfinder_run.log, whose lines read
#   Running: amrfinder -p bakta_output/<acc>/<acc>.faa --organism Acinetobacter_baumannii --plus --output amrfinder_results/<acc>.tsv --threads 2
# (protein mode on the Bakta .faa; the original loop script was not kept, so the loop below is a reconstruction).
# RGI v6.0.5 / CARD v4.0.1: RECONSTRUCTED from the output tables (protein mode, local DB, Perfect/Strict only; one .json +
# one .txt per genome, consistent with --clean). Alignment tool and threads were not recorded (RGI defaults used here).
# The batch loop itself is a reconstruction. Run from the project root.
set -euo pipefail
mkdir -p results/amr/amrfinder_results results/amr/card_results
while read -r ACC; do
  [[ -s results/amr/amrfinder_results/${ACC}.tsv ]] || \
    { echo "Running: amrfinder -p bakta_output/${ACC}/${ACC}.faa"; amrfinder -p bakta_output/${ACC}/${ACC}.faa --organism Acinetobacter_baumannii --plus \
              --output results/amr/amrfinder_results/${ACC}.tsv --threads 2; }
done < highquality_accessions.txt
# RGI (run in the rgi/card conda env)
while read -r ACC; do
  [[ -s results/amr/card_results/${ACC}.json ]] || \
    rgi main --input_sequence bakta_output/${ACC}/${ACC}.faa --input_type protein \
             --output_file results/amr/card_results/${ACC} --local --clean
done < highquality_accessions.txt
