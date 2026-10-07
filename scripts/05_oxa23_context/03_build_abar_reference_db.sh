#!/usr/bin/env bash
# 03_build_abar_reference_db.sh
# Builds the final (v4) curated BLAST database of characterised A. baumannii
# resistance islands / transposons used to classify blaOXA-23 genetic contexts.
#
# Why "isolated regions"? Early versions (v1-v3) used whole replicons
# (AB0057 chromosome NC_011586.2, Tn2008B source CP085788.1, Tn2009 source
# NZ_KM922672.1). Whole chromosomes/plasmids bias BLAST classification through
# broad background similarity rather than the transposon itself. Here the
# blaOXA-23-homologous region is located with the standalone Tn2006 probe
# (JN129846.1) and only that region (+5 kb flank each side) is kept.
#
# The isolation coordinates below were obtained with the locate step (Step 2)
# on the references downloaded from NCBI; re-run Step 2 if you change references.
set -euo pipefail

DB_DIR="databases/abar"
mkdir -p "${DB_DIR}"

echo ">> Step 1: Download reference sequences (NCBI nuccore, Entrez Direct)"
for acc in JX481978.1 JN129845.1 JN129847.1 HQ700358.2 JN129846.1 \
           EF059914.1 KP780408.1 CP085788.1 NZ_KM922672.1 NC_011586.2; do
    efetch -db nuccore -id "${acc}" -format fasta > "${DB_DIR}/${acc}.fna"
done
# AbaR25 JX481978.1 | TnAbaR4a JN129845.1 | AbaR4c fragment JN129847.1 | AbaR4 HQ700358.2
# Tn2006 JN129846.1 | Tn2007 EF059914.1 | Tn2008 KP780408.1
# Tn2008B CP085788.1 (whole chromosome) | Tn2009 NZ_KM922672.1 (whole plasmid)
# AB0057 NC_011586.2 (whole chromosome)

echo ">> Step 2: Locate the blaOXA-23 / Tn2006-homologous block in each whole replicon"
# Probe = standalone Tn2006 (JN129846.1); output sorted by bitscore.
for acc in CP085788.1 NZ_KM922672.1 NC_011586.2; do
    blastn -query "${DB_DIR}/JN129846.1.fna" -subject "${DB_DIR}/${acc}.fna" \
        -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
        2>/dev/null | sort -k10,10 -gr > "${DB_DIR}/hit_${acc}.tsv"
    echo "Top hits vs ${acc}:"; head -3 "${DB_DIR}/hit_${acc}.tsv"
done
# Observed top blocks used for this study:
#   CP085788.1    : 526,382-529,181  (~2.8 kb)  -> region 521,382-534,181 (block +/- 5 kb)
#   NZ_KM922672.1 :  69,129- 72,323  (~3.2 kb)  -> region  64,129- 77,323
#   NC_011586.2   : 582,113-589,181  (~7 kb)    -> region 577,113-594,181

echo ">> Step 3: Extract isolated regions (block + 5 kb flank each side)"
seqkit faidx "${DB_DIR}/CP085788.1.fna"    "CP085788.1:521382-534181"    > "${DB_DIR}/Tn2008B_isolated.fna"
seqkit faidx "${DB_DIR}/NZ_KM922672.1.fna" "NZ_KM922672.1:64129-77323"   > "${DB_DIR}/Tn2009_isolated.fna"
seqkit faidx "${DB_DIR}/NC_011586.2.fna"   "NC_011586.2:577113-594181"   > "${DB_DIR}/NC_011586_isolated.fna"
seqkit replace -p "^(.+)$" -r "CP085788.1_Tn2008B_isolated_region"   "${DB_DIR}/Tn2008B_isolated.fna"   -o "${DB_DIR}/Tn2008B_isolated_renamed.fna"
seqkit replace -p "^(.+)$" -r "NZ_KM922672.1_Tn2009_isolated_region" "${DB_DIR}/Tn2009_isolated.fna"    -o "${DB_DIR}/Tn2009_isolated_renamed.fna"
seqkit replace -p "^(.+)$" -r "NC_011586.2_AB0057_isolated_region"   "${DB_DIR}/NC_011586_isolated.fna" -o "${DB_DIR}/NC_011586_isolated_renamed.fna"
seqkit stats "${DB_DIR}"/*_isolated_renamed.fna     # expected ~12.8 kb, ~13.2 kb, ~17.1 kb

echo ">> Step 4: Assemble the final v4 database (10 sequences, no whole replicons)"
cat "${DB_DIR}/JX481978.1.fna" "${DB_DIR}/JN129845.1.fna" "${DB_DIR}/JN129847.1.fna" \
    "${DB_DIR}/HQ700358.2.fna" "${DB_DIR}/JN129846.1.fna" \
    "${DB_DIR}/NC_011586_isolated_renamed.fna" \
    "${DB_DIR}/EF059914.1.fna" "${DB_DIR}/KP780408.1.fna" \
    "${DB_DIR}/Tn2008B_isolated_renamed.fna" "${DB_DIR}/Tn2009_isolated_renamed.fna" \
    > "${DB_DIR}/AbaR_combined_v4.fna"
grep -c ">" "${DB_DIR}/AbaR_combined_v4.fna"        # expect 10
makeblastdb -in "${DB_DIR}/AbaR_combined_v4.fna" -dbtype nucl -out "${DB_DIR}/AbaR_db_v4"
seqkit fx2tab -nl "${DB_DIR}/AbaR_combined_v4.fna"  # reference lengths (used for normalised coverage)
