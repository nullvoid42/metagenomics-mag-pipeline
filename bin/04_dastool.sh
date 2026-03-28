#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════
# 04_dastool.sh — DAS Tool integration & dereplication
# ═══════════════════════════════════════════════════════════════════
# Usage: bash bin/04_dastool.sh <metabat2_dir> <maxbin2_dir> <concoct_dir> <outdir> <score_thresh> <samples>
set -euo pipefail

METABAT2_DIR="$1"; MAXBIN2_DIR="$2"; CONCOCT_DIR="$3"; OUTDIR="$4"; SCORE_THRESH="${5:-0.5}"; SAMPLES="$6"

mkdir -p "$OUTDIR/dereplicated_bins"

echo "═══════════════════════════════════════"
echo "  DAS Tool dereplication"
echo "  Score threshold: $SCORE_THRESH"
echo "═══════════════════════════════════════"

# ── Build FAB file lists ──────────────────────────────────────────
# DAS Tool needs a .tsv file listing binner outputs

# Build Scorerags for each binner
for SAMPLE_DIR in "$METABAT2_DIR"/*/; do
  [ -d "$SAMPLE_DIR" ] || continue
  SAMPLE=$(basename "$SAMPLE_DIR")
  mkdir -p "$OUTDIR/scores"
  
  # Create FAB file for metabat2
  BINNER_DIR="$METABAT2_DIR/$SAMPLE"
  FFILE="$OUTDIR/scores/metabat2_${SAMPLE}.txt"
  echo -e "bin_id\tfile" > "$FFILE"
  for bin in "$BINNER_DIR"/${SAMPLE}.bin.*.fa; do
    [ -f "$bin" ] || continue
    echo -e "$(basename $bin)\t$bin" >> "$FFILE"
  done
  
  # Create FAB file for maxbin2
  BINNER_DIR="$MAXBIN2_DIR/$SAMPLE"
  FFILE="$OUTDIR/scores/maxbin2_${SAMPLE}.txt"
  echo -e "bin_id\tfile" > "$FFILE"
  for bin in "$BINNER_DIR"/${SAMPLE}.maxbin.*.fasta; do
    [ -f "$bin" ] || continue
    echo -e "$(basename $bin)\t$bin" >> "$FFILE"
  done
  
  # Create FAB file for concoct
  BINNER_DIR="$CONCOCT_DIR/$SAMPLE"
  FFILE="$OUTDIR/scores/concoct_${SAMPLE}.txt"
  echo -e "bin_id\tfile" > "$FFILE"
  for bin in "$BINNER_DIR"/bin_*.fa; do
    [ -f "$bin" ] || continue
    echo -e "$(basename $bin)\t$bin" >> "$FFILE"
  done
done

# ── Run DAS_Tool ──────────────────────────────────────────────────
# Get contigs from first sample's assembly
SAMPLE_FIRST=$(echo "$SAMPLES" | awk '{print $1}')
CONTIGS="results/02_assembly/$SAMPLE_FIRST/${SAMPLE_FIRST}.contigs.fa"

echo "  Running DAS Tool..."
DAS_TOOL_CMD="Fasta_To_Scaffolding_Collection_Generator.pl \
  $CONTIGS \
  $OUTDIR/scores/metabat2_*.txt \
  $OUTDIR/scores/maxbin2_*.txt \
  $OUTDIR/scores/concoct_*.txt \
  -o $OUTDIR/scores/scaffolds2bin.tsv"

eval "$DAS_TOOL_CMD"

# Run DAS_Tool
DAS_Tool \
  --bins "$OUTDIR/scores/scaffolds2bin.tsv" \
  --output_dir "$OUTDIR" \
  --search_engine blast \
  --score_threshold "$SCORE_THRESH" \
  --create_plots \
  --threads 16

# Move dereplicated bins
if [ -d "$OUTDIR/DASTool_bins" ]; then
  mv "$OUTDIR/DASTool_bins"/* "$OUTDIR/dereplicated_bins/"
  rm -rf "$OUTDIR/DASTool_bins"
fi

N_BINS=$(ls "$OUTDIR/dereplicated_bins"/*.fa 2>/dev/null | wc -l | tr -d ' ')
echo "─── DAS Tool Summary ──────────────────"
echo "  Dereplicated bins: $N_BINS"
echo "═══════════════════════════════════════"
echo "✓ DAS Tool complete: $N_BINS high-quality bins"
