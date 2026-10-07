#!/bin/bash
# Verified script: copied unchanged from the cluster (scripts_logs/run_bakta.sh).
# Bakta v1.11.3, database v6.0 (full). Run from the project root.
# Annotates every genome in highquality_fastas/ (7 parallel jobs x 8 threads, up to 3 attempts each,
# skips genomes that already have a .json). Note: the folder held 898 FASTAs when it was run (897 study
# genomes + GCF_009035845.1, which was later excluded; Panaroo used the 897 in highquality_accessions.txt).
mkdir -p bakta_output

run_one() {
    fasta="$1"
    sample=$(basename "$fasta" .fna)
    if [[ -f "bakta_output/${sample}/${sample}.json" ]]; then
        echo "Skipping ${sample}, already done"
        return
    fi

    for attempt in 1 2 3; do
        bakta --db /home/zbecerra/bakta_db/db \
              --output bakta_output/${sample} \
              --prefix ${sample} \
              --threads 8 \
              --force \
              "$fasta" && break
        echo "Attempt $attempt failed for ${sample}, retrying..."
    done
}
export -f run_one

ls highquality_fastas/*.fna | xargs -P 7 -I {} bash -c 'run_one "$@"' _ {}
