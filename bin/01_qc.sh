#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════
# 01_qc.sh — Quality control with fastp
# ═══════════════════════════════════════════════════════════════════
# Usage: bash bin/01_qc.sh <R1.fastq.gz> <R2.fastq.gz> <outdir> <sample> [adapters]
set -euo pipefail

R1="$1"; R2="$2"; OUTDIR="$3"; SAMPLE="$4"; ADAPTERS="${5:-}"

# ── Validate inputs ──────────────────────────────────────────────
if [ ! -f "$R1" ]; then echo "ERROR: R1 not found: $R1"; exit 1; fi
if [ ! -f "$R2" ]; then echo "ERROR: R2 not found: $R2"; exit 1; fi

mkdir -p "$OUTDIR" "$OUTDIR/unpaired"

echo "═══════════════════════════════════════"
echo "  fastp QC for: $SAMPLE"
echo "  R1: $R1"
echo "  R2: $R2"
echo "═══════════════════════════════════════"

# ── Build fastp command ───────────────────────────────────────────
CMD="fastp \
  --in1 '$R1' --in2 '$R2' \
  --out1 '$OUTDIR/${SAMPLE}_1.filtered.fastq.gz' \
  --out2 '$OUTDIR/${SAMPLE}_2.filtered.fastq.gz' \
  --json '$OUTDIR/${SAMPLE}.fastp.json' \
  --html '$OUTDIR/${SAMPLE}.fastp.html' \
  --detect_adapter_for_pe \
  --cut_front --cut_tail \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --length_required 50 \
  --low_complexity_filter \
  --thread 8 \
  --report_title '$SAMPLE fastp report'"

# Add custom adapters if provided
if [ -n "$ADAPTERS" ]; then
  CMD="$CMD --adapter_fasta '$ADAPTERS'"
fi

# ── Run ──────────────────────────────────────────────────────────
eval "$CMD"

echo "✓ QC complete: $OUTDIR/${SAMPLE}_1.filtered.fastq.gz"
echo "✓ Report: $OUTDIR/${SAMPLE}.fastp.html"
echo "✓ JSON: $OUTDIR/${SAMPLE}.fastp.json"

# ── Quick stats ───────────────────────────────────────────────────
READS_RAW=$(zcat "$R1" | wc -l | awk '{print $1/4}')
READS_FILT=$(zcat "$OUTDIR/${SAMPLE}_1.filtered.fastq.gz" | wc -l | awk '{print $1/4}')
echo "─── QC Summary ───────────────────────"
echo "  Raw reads:    $READS_RAW"
echo "  Filtered:     $READS_FILT"
echo "  Retained:     $(echo "scale=2; $READS_FILT/$READS_RAW*100" | bc)%"
echo "═══════════════════════════════════════"
