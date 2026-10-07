#!/bin/bash
# Idempotent patch for mlst 2.17.6 (Seemann): the script compares the BLAST+ version as a decimal, so
# "2.16 >= 2.9" is false in Perl and mlst stops with "needs BLAST+ 2.9.0 or later".
# This removes that check. Run it inside the conda environment that holds mlst, and again after cloning/renaming the environment
# (the patch lives in <env>/bin/mlst and is not carried over). Check first whether a newer mlst release already fixes this.
# Usage: bash scripts/00_setup/patch_mlst_blast_check.sh [path/to/mlst]
set -euo pipefail
MLST="${1:-$(which mlst)}"
if grep -q 'decver >= 2\.9 or err' "$MLST"; then
  [ -e "$MLST.bak" ] || cp "$MLST" "$MLST.bak"
  sed -i 's/\$decver >= 2\.9 or err/1 or err/' "$MLST"
  echo "patched: $MLST (backup: $MLST.bak)"
else
  echo "nothing to do: $MLST has no BLAST+ version check to patch (already patched or fixed upstream)"
fi
