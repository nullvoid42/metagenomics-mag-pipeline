#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════
# 03c_concoct.sh — CONCOCT binning
# ═══════════════════════════════════════════════════════════════════
# Usage: bash bin/03c_concoct.sh <contigs.fa> <outdir> <min_contig>
set -euo pipefail

CONTIGS="$1"; OUTDIR="$2"; MIN_CONTIG="${3:-1500}"

SAMPLE=$(basename "$(dirname "$OUTDIR")")
CONTIGS_BASE=$(basename "$CONTIGS" .fa)

mkdir -p "$OUTDIR"

echo "═══════════════════════════════════════"
echo "  CONCOCT binning: $SAMPLE"
echo "  Contigs: $CONTIGS"
echo "═══════════════════════════════════════"

QC_DIR="results/01_qc"
R1="${QC_DIR}/${SAMPLE}_1.filtered.fastq.gz"
R2="${QC_DIR}/${SAMPLE}_2.filtered.fastq.gz"

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
  echo "ERROR: Filtered reads not found. Run QC first."
  exit 1
fi

CONTIGS_BASE_NAME=$(basename "$CONTIGS" .fa)

# ── Step 1: cut_up_fasta.py ──────────────────────────────────────
echo "  Cutting contigs..."
cut_up_fasta.py "$CONTIGS" -c 10000 -o 0 | sed 's/ />custom_/g' > "$OUTDIR/contigs_10K.fasta"

# ── Step 2: bowtie2 mapping ──────────────────────────────────────
ASSEMBLY_DIR=$(dirname "$CONTIGS")
BT2_IDX="$ASSEMBLY_DIR/${CONTIGS_BASE}.bt2"

# Build index if not exists
if [ ! -f "${BT2_IDX}.1.bt2" ]; then
  echo "  Building bowtie2 index..."
  bowtie2-build "$CONTIGS" "$BT2_IDX" -p 16
fi

echo "  Mapping reads..."
BEDFILE="$OUTDIR/${SAMPLE}.bed"
BAMFILE="$OUTDIR/${SAMPLE}.sorted.bam"
bowtie2 -x "$BT2_IDX" -1 "$R1" -2 "$R2" -p 16 --no-unal \
  | samtools sort -@ 8 -o "$BAMFILE" -
bedtools genomecov -ibam "$BAMFILE" -bga \
  | awk '$4>0' > "$BEDFILE"
rm -f "$BAMFILE"

# ── Step 3: concoct_coverage_table.py ────────────────────────────
echo "  Computing coverage table..."
concoct_coverage_table.py "$BEDFILE" "$OUTDIR/contigs_10K.fasta" > "$OUTDIR/coverage_table.tsv"

# ── Step 4: concoct ────────────────────────────────────────────────
echo "  Running CONCOCT..."
concoct \
  --composition_file "$OUTDIR/contigs_10K.fasta" \
  --coverage_table "$OUTDIR/coverage_table.tsv" \
  --output_directory "$OUTDIR" \
  -t 16 \
  -c 40 \
  -l "${MIN_CONTIG}"

# ── Step 5: merge_cutup_clustering.py ───────────────────────────
echo "  Merging clusterings..."
merge_cutup_clustering.py "$OUTDIR/clustering_gt${MIN_CONTIG}.csv" > "$OUTDIR/clustering_merged.csv"

# ── Step 6: extract_fasta_bins.py ────────────────────────────────
echo "  Extracting bins..."
extract_fasta_bins.py "$CONTIGS" "$OUTDIR/clustering_merged.csv" \
  --output_path "$OUTDIR" \
  --gzip_format

N_BINS=$(ls "$OUTDIR"/bin_*.fa 2>/dev/null | wc -l | tr -d ' ')
echo "─── CONCOCT Summary ───────────────────"
echo "  Bins recovered: $N_BINS"
echo "═══════════════════════════════════════"
echo "✓ CONCOCT complete: $N_BINS bins"
