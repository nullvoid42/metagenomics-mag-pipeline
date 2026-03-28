#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════
# 02_assemble.sh — MEGAHIT assembly
# ═══════════════════════════════════════════════════════════════════
# Usage: bash bin/02_assemble.sh <R1.filtered.fastq.gz> <R2.filtered.fastq.gz> <outdir> <sample> <min_contig>
set -euo pipefail

R1="$1"; R2="$2"; OUTDIR="$3"; SAMPLE="$4"; MIN_CONTIG="${5:-1500}"

mkdir -p "$OUTDIR"

echo "═══════════════════════════════════════"
echo "  MEGAHIT assembly: $SAMPLE"
echo "  R1: $R1"
echo "  R2: $R2"
echo "  Outdir: $OUTDIR"
echo "═══════════════════════════════════════"

# Check for proxies (common on HPC)
PROXY_ENV=""
if [ -n "${http_proxy:-}" ]; then
  PROXY_ENV="http_proxy=$http_proxy https_proxy=$https_proxy"
fi

# MEGAHIT parameters tuned for metagenomes
$PROXY_ENV megahit \
  -1 "$R1" \
  -2 "$R2" \
  -o "$OUTDIR" \
  -t 32 \
  --min-count 2 \
  --k-min 27 \
  --k-max 91 \
  --k-step 10 \
  --memory 0.4 \
  --presets meta-large \
  2>&1 | tee "$OUTDIR/megahit.log"

# Rename final contigs file
if [ -f "$OUTDIR/final.contigs.fa" ]; then
  mv "$OUTDIR/final.contigs.fa" "$OUTDIR/${SAMPLE}.contigs.fa"
  echo "✓ Renamed final.contigs.fa → ${SAMPLE}.contigs.fa"
fi

CONTIGS="$OUTDIR/${SAMPLE}.contigs.fa"
if [ ! -f "$CONTIGS" ]; then
  echo "ERROR: Contigs file not created: $CONTIGS"
  exit 1
fi

# Filter contigs by minimum length
CONTIGS_FILTERED="$OUTDIR/${SAMPLE}.contigs.min${MIN_CONTIG}.fa"
awk -v min="$MIN_CONTIG" '/^>/{header=$0; getline seq; if(length(seq)>=min){print header"\n"seq}}' "$CONTIGS" > "$CONTIGS_FILTERED"

N_CONTIGS=$(grep -c "^>" "$CONTIGS_FILTERED")
TOTAL_BP=$(awk '/^[^>]/{sum+=length($0)}END{print sum}' "$CONTIGS_FILTERED")
echo "─── Assembly Summary ──────────────────"
echo "  Sample:       $SAMPLE"
echo "  Contigs ≥ ${MIN_CONTIG}bp: $N_CONTIGS"
echo "  Total bp:     $TOTAL_BP"
echo "  Avg length:   $(echo "scale=0; $TOTAL_BP/$N_CONTIGS" | bc 2>/dev/null || echo 'N/A')"
echo "  N50:          $(awk -v min="$MIN_CONTIG" '/^>/{getline seq; if(length(seq)>=min) print length(seq)}' "$CONTIGS_FILTERED" | sort -n | awk '{a[NR]=$1; s+=$1} END{print a[int(NR*0.5)]}')"
echo "═══════════════════════════════════════"
echo "✓ Assembly complete: $CONTIGS_FILTERED"
