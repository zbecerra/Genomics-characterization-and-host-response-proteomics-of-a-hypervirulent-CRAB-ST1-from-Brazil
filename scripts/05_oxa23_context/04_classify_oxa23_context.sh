#!/usr/bin/env bash
# 04_classify_oxa23_context.sh
# Classifies every individually extracted blaOXA-23 flanking region (all 848
# copies, not only cluster representatives) against AbaR_db_v4.
#
# Per query: single best hit (highest query coverage, then bitscore).
# Coverage is normalised to the SHORTER of query/reference so that complete
# matches to short references (e.g. Tn2008, 4,593 bp) are not penalised:
#       normalised coverage = alignment_length / min(qlen, slen) x 100
# Queries with < 99 % identity OR < 90 % normalised coverage -> "Unassigned/atypical".
# AbaR4, Tn2006 and the isolated AB0057 region overlap in sequence and are
# collapsed into a single "Tn2006/AbaR4-type" category.
#
# Result obtained on the 848 copies of this study:
#   Tn2006/AbaR4-type 538 | Unassigned/atypical 242 | Tn2009 37 | Tn2008 26 |
#   Tn2007 3 | AbaR25 2 | Tn2008B 0
set -euo pipefail

DB="databases/abar/AbaR_db_v4"
QUERY="oxa23_context/oxa23_contexts.fasta"
OUT="oxa23_context"
MIN_ID=99
MIN_NCOV=90

blastn -query "${QUERY}" -db "${DB}" \
       -outfmt "6 qseqid sseqid pident length qcovs qlen slen evalue bitscore" \
       -max_target_seqs 1 -out "${OUT}/all_copies_vs_AbaR_db_v4_full.tsv"

# best hit per query: qcovs (col5) then bitscore (col9)
sort -k1,1 -k5,5 -gr -k9,9 -gr "${OUT}/all_copies_vs_AbaR_db_v4_full.tsv" | \
    awk '!seen[$1]++' > "${OUT}/best_hit_v4_full.tsv"

# normalised coverage table: query, subject, pident, norm_cov
awk -F'\t' '{minlen=($6<$7)?$6:$7; printf "%s\t%s\t%.2f\t%.2f\n",$1,$2,$3,($4/minlen)*100}' \
    "${OUT}/best_hit_v4_full.tsv" > "${OUT}/normalized_coverage.tsv"

# final per-copy classification
awk -F'\t' -v mi="${MIN_ID}" -v mc="${MIN_NCOV}" 'BEGIN{OFS="\t"; print "copy","best_reference","pident","norm_cov","class"}
{
  cls="Unassigned/atypical"
  if ($3>=mi && $4>=mc) {
    if      ($2 ~ /^HQ700358.2|^JN129846.1|^NC_011586.2/) cls="Tn2006/AbaR4-type"
    else if ($2 ~ /^JN129845.1/)                          cls="TnAbaR4a"
    else if ($2 ~ /^JN129847.1/)                          cls="AbaR4c"
    else if ($2 ~ /^JX481978.1/)                          cls="AbaR25"
    else if ($2 ~ /^EF059914.1/)                          cls="Tn2007"
    else if ($2 ~ /^KP780408.1/)                          cls="Tn2008"
    else if ($2 ~ /^CP085788.1/)                          cls="Tn2008B"
    else if ($2 ~ /^NZ_KM922672.1/)                       cls="Tn2009"
  }
  print $1,$2,$3,$4,cls
}' "${OUT}/normalized_coverage.tsv" > "${OUT}/oxa23_context_classification.tsv"

echo "=== Classification summary (copies) ==="
awk -F'\t' 'NR>1{print $5}' "${OUT}/oxa23_context_classification.tsv" | sort | uniq -c | sort -rn
