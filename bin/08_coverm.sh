#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════
# 08_coverm.sh — CoverM MAG abundance estimation
# ═══════════════════════════════════════════════════════════════════
# Usage: bash bin/08_coverm.sh <bins_dir> <outdir> <threads> <samples>
set -euo pipefail

BINS_DIR="$1"; OUTDIR="$2"; THREADS="${3:-16}"; SAMPLES="$4"

mkdir -p "$OUTDIR"

echo "═══════════════════════════════════════"
echo "  CoverM abundance estimation"
echo "  Bins dir: $BINS_DIR"
echo "  Samples: $SAMPLES"
echo "═══════════════════════════════════════"

N_BINS=$(find "$BINS_DIR" -name "*.fa" 2>/dev/null | wc -l | tr -d ' ')
echo "  Reference MAGs: $N_BINS"

# ── Build reference database from bins ───────────────────────────
REF_FASTA="$OUTDIR/reference_mags.fa"
cat "$BINS_DIR"/*.fa > "$REF_FASTA" 2>/dev/null || true

# Index with CoverM's built-in method (or use minimap2)
echo "  Building CoverM reference..."

# ── Run CoverM on each sample ────────────────────────────────────
COVERM_OUT="$OUTDIR/coverm_raw.tsv"
SAMPLE_ARRAY=($SAMPLES)

R1_FILES=""
R2_FILES=""
for SAMPLE in "${SAMPLE_ARRAY[@]}"; do
  R1_FILES="$R1_FILES results/01_qc/${SAMPLE}_1.filtered.fastq.gz"
  R2_FILES="$R2_FILES results/01_qc/${SAMPLE}_2.filtered.fastq.gz"
done

# CoverM genome mode - compute coverage directly
coverm genome \
  --genome-fasta-directory "$BINS_DIR" \
  --reference "$REF_FASTA" \
  --interleaved $(echo $R1_FILES | awk '{print $1}') \
  --threads "$THREADS" \
  --min-read-percent-identity 0.95 \
  --min-read-alignment-length 45 \
  --trim-upper 0.0 \
  --trim-lower 0.1 \
  --output-file "$COVERM_OUT" \
  --output-format tsv \
  2>&1 | tee "$OUTDIR/coverm.log"

# ── Generate summary with taxonomy ────────────────────────────────
GTDBTK_SUMMARY="results/06_gtdbtk/gtdbtk_summary.tsv"
SUMMARY_OUT="$OUTDIR/coverm_summary.tsv"

if [ -f "$GTDBTK_SUMMARY" ]; then
  python3 -c "
import sys
# Read taxonomy
taxonomy = {}
with open('$GTDBTK_SUMMARY') as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 8: taxonomy[parts[0]] = parts[7]  # genus

# Read coverm
print('mag\tgenus\tabundance\tcoverage\tlength\tcount')
with open('$COVERM_OUT') as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 5: continue
        mag = parts[0].replace('.fa','')
        genus = taxonomy.get(mag, 'Unknown')
        abundance = parts[1]
        coverage = parts[2]
        length = parts[3]
        count = parts[4]
        print(f'{mag}\t{genus}\t{abundance}\t{coverage}\t{length}\t{count}')
" > "$SUMMARY_OUT"
else
  cp "$COVERM_OUT" "$SUMMARY_OUT"
fi

echo "─── CoverM Summary ───────────────────"
echo "  MAGs: $N_BINS"
echo "  Output: $SUMMARY_OUT"
echo "═══════════════════════════════════════"
echo "✓ CoverM complete"
