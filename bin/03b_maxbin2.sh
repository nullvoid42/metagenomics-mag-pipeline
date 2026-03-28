#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════
# 03b_maxbin2.sh — MaxBin2 binning
# ═══════════════════════════════════════════════════════════════════
# Usage: bash bin/03b_maxbin2.sh <contigs.fa> <outdir> <min_contig>
set -euo pipefail

CONTIGS="$1"; OUTDIR="$2"; MIN_CONTIG="${3:-1500}"

SAMPLE=$(basename "$(dirname "$OUTDIR")")
CONTIGS_BASE=$(basename "$CONTIGS" .fa)

mkdir -p "$OUTDIR"

echo "═══════════════════════════════════════"
echo "  MaxBin2 binning: $SAMPLE"
echo "  Contigs: $CONTIGS"
echo "═══════════════════════════════════════"

# Get filtered reads
QC_DIR="results/01_qc"
R1="${QC_DIR}/${SAMPLE}_1.filtered.fastq.gz"
R2="${QC_DIR}/${SAMPLE}_2.filtered.fastq.gz"

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
  echo "ERROR: Filtered reads not found. Run QC first."
  exit 1
fi

# Run MaxBin2
run_MaxBin.pl \
  -thread 16 \
  -contig "$CONTIGS" \
  -reads1 "$R1" \
  -reads2 "$R2" \
  -out "$OUTDIR/${SAMPLE}.maxbin" \
  -min_contig_length "${MIN_CONTIG}"

N_BINS=$(ls "$OUTDIR"/${SAMPLE}.maxbin.*.fasta 2>/dev/null | wc -l | tr -d ' ')
echo "─── MaxBin2 Summary ───────────────────"
echo "  Bins recovered: $N_BINS"
echo "═══════════════════════════════════════"
echo "✓ MaxBin2 complete: $N_BINS bins"
