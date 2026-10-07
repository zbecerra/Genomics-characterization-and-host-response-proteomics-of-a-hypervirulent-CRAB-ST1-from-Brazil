#!/bin/bash
# =============================================================================
# Pipeline: Download Acinetobacter baumannii ST1 genomes ONE BY ONE
#   - Prefer GCF (RefSeq) over GCA (GenBank)
#   - Keep only latest version
#   - Check MLST (Pasteur scheme, abaumannii_2) → keep ST1 only
#   - Survives terminal disconnect via nohup
# =============================================================================

set -euo pipefail

OUTDIR="ACB_ST1_pipeline"
TARGET_ST="1"
SCHEME="abaumannii_2"   # Pasteur scheme (corrected; see 02_reverify_mlst_pasteur.sh)
TAXID="470"  # Acinetobacter baumannii

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; CYAN='\033[0;36m'; NC='\033[0m'
info()  { echo -e "${GREEN}[INFO]${NC}    $*"; }
warn()  { echo -e "${YELLOW}[REMOVED]${NC} $*"; }
kept()  { echo -e "${CYAN}[KEPT]${NC}    $*"; }
error() { echo -e "${RED}[ERROR]${NC}   $*"; exit 1; }

for cmd in datasets mlst; do
  command -v "$cmd" &>/dev/null || error "'$cmd' not found"
done

mkdir -p "$OUTDIR"/{accessions,kept_ST1,temp,logs}
cd "$OUTDIR"

PROCESSED_LOG="logs/processed.txt"
MLST_LOG="logs/mlst_results_all.tsv"
touch "$PROCESSED_LOG"
[[ ! -f "$MLST_LOG" ]] && echo -e "accession\tscheme\tST\talleles" > "$MLST_LOG"

# =============================================================================
# STEP 1 — Fetch dehydrated manifest
# =============================================================================
if [[ ! -f accessions/accessions.txt ]]; then
  info "Step 1: Downloading dehydrated manifest (Acinetobacter baumannii, taxid $TAXID)..."
  
  datasets download genome taxon "$TAXID" \
    --dehydrated \
    --filename accessions/acb_dehydrated.zip \
    --no-progressbar

  info "Unzipping..."
  unzip -o accessions/acb_dehydrated.zip -d accessions/dehydrated/ > /dev/null 2>&1

  info "Extracting accessions..."
  grep -o 'GC[AF]_[0-9]*\.[0-9]*' accessions/dehydrated/ncbi_dataset/fetch.txt \
    | sort -u > accessions/all_accessions.txt

  TOTAL_FOUND=$(wc -l < accessions/all_accessions.txt)
  info "Found $TOTAL_FOUND accessions. Deduplicating..."
  
  python3 << 'PYTHON'
import re
accessions = {}
with open('accessions/all_accessions.txt') as f:
    for line in f:
        acc = line.strip()
        if not acc: continue
        m = re.match(r'(GC[AF])_(\d+)\.(\d+)', acc)
        if not m: continue
        db_type, acc_id, version = m.groups()
        version = int(version)
        key = (acc_id, version)
        if key not in accessions:
            accessions[key] = (db_type, acc)
        else:
            existing_type, existing_acc = accessions[key]
            if db_type == 'GCF' and existing_type == 'GCA':
                accessions[key] = (db_type, acc)

final = {}
for (acc_id, version), (db_type, acc) in accessions.items():
    if acc_id not in final:
        final[acc_id] = (version, acc)
    else:
        existing_version, existing_acc = final[acc_id]
        if version > existing_version:
            final[acc_id] = (version, acc)

with open('accessions/accessions.txt', 'w') as f:
    for acc_id in sorted(final.keys()):
        _, acc = final[acc_id]
        f.write(acc + '\n')
print(f"Deduplicated: {len(final)} unique genomes")
PYTHON

  info "Ready to process $(wc -l < accessions/accessions.txt) genomes"
else
  info "Accession list exists ($(wc -l < accessions/accessions.txt) genomes). Skipping fetch."
fi

TOTAL=$(wc -l < accessions/accessions.txt)
echo ""
info "Step 2: Processing $TOTAL genomes..."
echo ""

KEPT_COUNT=0
REMOVED_COUNT=0
FAILED_COUNT=0
IDX=0

while read -r ACCESSION; do
  [[ -z "$ACCESSION" ]] && continue
  IDX=$((IDX + 1))

  if grep -qx "$ACCESSION" "$PROCESSED_LOG"; then
    info "[$IDX/$TOTAL] $ACCESSION — already done"
    continue
  fi

  echo "── [$IDX/$TOTAL] $ACCESSION"
  TMPDIR="temp/${ACCESSION}"
  mkdir -p "$TMPDIR"

  if ! datasets download genome accession "$ACCESSION" \
      --include genome \
      --filename "${TMPDIR}/genome.zip" \
      --no-progressbar 2>>"logs/download_errors.txt"; then
    warn "[$IDX/$TOTAL] $ACCESSION — download failed"
    FAILED_COUNT=$((FAILED_COUNT + 1))
    echo "$ACCESSION" >> "$PROCESSED_LOG"
    rm -rf "$TMPDIR"
    continue
  fi

  unzip -o "${TMPDIR}/genome.zip" -d "${TMPDIR}/raw" > /dev/null 2>&1
  FASTA=$(find "${TMPDIR}/raw" -name "*.fna" | head -1)

  if [[ -z "$FASTA" ]]; then
    warn "[$IDX/$TOTAL] $ACCESSION — no FASTA"
    FAILED_COUNT=$((FAILED_COUNT + 1))
    echo "$ACCESSION" >> "$PROCESSED_LOG"
    rm -rf "$TMPDIR"
    continue
  fi

  MLST_OUT=$(mlst --scheme "$SCHEME" --quiet "$FASTA" 2>/dev/null || true)
  ST=$(echo "$MLST_OUT" | awk -F'\t' '{print $3}')
  echo -e "${ACCESSION}\t${MLST_OUT#*$'\t'}" >> "$MLST_LOG"

  if [[ "$ST" == "$TARGET_ST" ]]; then
    cp "$FASTA" "kept_ST1/${ACCESSION}.fna"
    kept "[$IDX/$TOTAL] $ACCESSION — ST${ST} ✓"
    KEPT_COUNT=$((KEPT_COUNT + 1))
  else
    warn "[$IDX/$TOTAL] $ACCESSION — ST${ST}"
    REMOVED_COUNT=$((REMOVED_COUNT + 1))
  fi

  rm -rf "$TMPDIR"
  echo "$ACCESSION" >> "$PROCESSED_LOG"

  if (( IDX % 50 == 0 )); then
    echo ""
    info "Progress: $IDX/$TOTAL | Kept: $KEPT_COUNT | Removed: $REMOVED_COUNT | Failed: $FAILED_COUNT"
    echo ""
  fi

done < accessions/accessions.txt

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  COMPLETE"
echo "════════════════════════════════════════════════════════════════"
echo "  Processed: $IDX | Kept: $KEPT_COUNT | Removed: $REMOVED_COUNT | Failed: $FAILED_COUNT"
echo "  ST1 genomes → $OUTDIR/kept_ST1/"
echo "════════════════════════════════════════════════════════════════"
