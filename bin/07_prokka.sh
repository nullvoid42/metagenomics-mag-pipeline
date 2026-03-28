#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════
# 07_prokka.sh — Prokka gene prediction & annotation
# ═══════════════════════════════════════════════════════════════════
# Usage: bash bin/07_prokka.sh <bins_dir> <outdir> <kingdom> <threads>
set -euo pipefail

BINS_DIR="$1"; OUTDIR="$2"; KINGDOM="${3:-Bacteria}"; THREADS="${4:-16}"

mkdir -p "$OUTDIR"

echo "═══════════════════════════════════════"
echo "  Prokka gene annotation"
echo "  Bins dir: $BINS_DIR"
echo "  Kingdom: $KINGDOM"
echo "═══════════════════════════════════════"

N_BINS=$(find "$BINS_DIR" -name "*.fa" -o -name "*.fasta" 2>/dev/null | wc -l | tr -d ' ')
echo "  Annotating $N_BINS MAGs..."

COUNTER=0
for bin in "$BINS_DIR"/*.fa; do
  [ -f "$bin" ] || continue
  BIN_NAME=$(basename "$bin" .fa)
  COUNTER=$((COUNTER+1))
  
  echo "  [$COUNTER/$N_BINS] Annotating: $BIN_NAME"
  
  prokka \
    "$bin" \
    --outdir "$OUTDIR/$BIN_NAME" \
    --prefix "$BIN_NAME" \
    --kingdom "$KINGDOM" \
    --cpus "$THREADS" \
    --add_genes \
    --compliant \
    --force \
    2>&1 | grep -v "^$"
  
  # Copy to combined directory
  mkdir -p "$OUTDIR/all_annotations"
  cp "$OUTDIR/$BIN_NAME/$BIN_NAME.gbk" "$OUTDIR/all_annotations/" 2>/dev/null || true
  cp "$OUTDIR/$BIN_NAME/$BIN_NAME.gff" "$OUTDIR/all_annotations/" 2>/dev/null || true
done

echo "─── Prokka Summary ───────────────────"
echo "  Annotated: $COUNTER MAGs"
echo "  Output: $OUTDIR/"
echo "═══════════════════════════════════════"
echo "✓ Prokka complete"
