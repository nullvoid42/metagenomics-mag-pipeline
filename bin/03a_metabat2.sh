#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════
# 03a_metabat2.sh — MetaBAT2 binning
# ═══════════════════════════════════════════════════════════════════
# Usage: bash bin/03a_metabat2.sh <contigs.fa> <outdir> <min_contig>
set -euo pipefail

CONTIGS="$1"; OUTDIR="$2"; MIN_CONTIG="${3:-1500}"

SAMPLE=$(basename "$(dirname "$OUTDIR")")
CONTIGS_BASE=$(basename "$CONTIGS" .fa)
CONTIGS_FILTERED="${CONTIGS%.fa}.min${MIN_CONTIG}.fa"

mkdir -p "$OUTDIR"

echo "═══════════════════════════════════════"
echo "  MetaBAT2 binning: $SAMPLE"
echo "  Contigs: $CONTIGS_FILTERED"
echo "═══════════════════════════════════════"

# ── Compute coverage with bowtie2 + jgi_summarize_bam_contig_depths ─────────
CONTIGS_SHORT="${CONTIGS_FILTERED%.fa}.min${MIN_CONTIG}.fa"

# Get read files from QC output
QC_DIR="results/01_qc"
R1="${QC_DIR}/${SAMPLE}_1.filtered.fastq.gz"
R2="${QC_DIR}/${SAMPLE}_2.filtered.fastq.gz"

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
  echo "ERROR: Filtered reads not found. Run QC first."
  echo "  Expected: $R1, $R2"
  exit 1
fi

# Build bowtie2 index (one-time, cache in assembly dir)
ASSEMBLY_DIR=$(dirname "$CONTIGS")
BT2_IDX="$ASSEMBLY_DIR/${CONTIGS_BASE}.bt2"
if [ ! -f "${BT2_IDX}.1.bt2" ]; then
  echo "  Building bowtie2 index..."
  bowtie2-build "$CONTIGS" "$BT2_IDX" -p 16
fi

# Map reads
BAM="$OUTDIR/${SAMPLE}.sorted.bam"
echo "  Mapping reads to contigs..."
bowtie2 -x "$BT2_IDX" -1 "$R1" -2 "$R2" -p 16 --no-unal \
  | samtools sort -@ 8 -o "$BAM" -

# ── jgi_summarize_depths ──────────────────────────────────────────
DEPTH_FILE="$OUTDIR/${SAMPLE}.depth.txt"
jgi_summarize_bam_contig_depths --outputDepth "$DEPTH_FILE" "$BAM"

# Remove temp BAM to save space
rm -f "$BAM"

# ── MetaBAT2 binning ──────────────────────────────────────────────
echo "  Running MetaBAT2..."
metabat2 \
  -i "$CONTIGS" \
  -a "$DEPTH_FILE" \
  -o "$OUTDIR/${SAMPLE}.bin" \
  -t 16 \
  --minContig "${MIN_CONTIG}" \
  --minCV 0.1 \
  --minCVSum 0.1 \
  --maxP 95 \
  --minS 0.6 \
  --maxEdges 200

N_BINS=$(ls "$OUTDIR"/${SAMPLE}.bin.*.fa 2>/dev/null | wc -l | tr -d ' ')
echo "─── MetaBAT2 Summary ───────────────────"
echo "  Bins recovered: $N_BINS"
echo "═══════════════════════════════════════"
echo "✓ MetaBAT2 complete: $N_BINS bins"
