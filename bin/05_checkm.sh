#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════
# 05_checkm.sh — CheckM quality assessment
# ═══════════════════════════════════════════════════════════════════
# Usage: bash bin/05_checkm.sh <bins_dir> <outdir> <lineage_wf> <threads>
set -euo pipefail

BINS_DIR="$1"; OUTDIR="$2"; LINEAGE_WF="${3:-true}"; THREADS="${4:-16}"

mkdir -p "$OUTDIR"

echo "═══════════════════════════════════════"
echo "  CheckM quality assessment"
echo "  Bins dir: $BINS_DIR"
echo "  Lineage workflow: $LINEAGE_WF"
echo "═══════════════════════════════════════"

# Check if bins exist
N_BINS=$(find "$BINS_DIR" -name "*.fa" -o -name "*.fasta" 2>/dev/null | wc -l | tr -d ' ')
if [ "$N_BINS" -eq 0 ]; then
  echo "ERROR: No bins found in $BINS_DIR"
  exit 1
fi
echo "  Found $N_BINS bins"

# Run CheckM
if [ "$LINEAGE_WF" = "true" ]; then
  # Lineage-specific markers (more accurate)
  checkm lineage_wf \
    "$BINS_DIR" \
    "$OUTDIR/checkm_lineage" \
    -x fa \
    -t "$THREADS" \
    --reduced_tree \
    2>&1 | tee "$OUTDIR/checkm.log"
  
  # Generate summary
  checkm qa \
    "$OUTDIR/checkm_lineage/checkm.lineage.ms" \
    "$OUTDIR/checkm_lineage" \
    -o 2 \
    -f "$OUTDIR/checkm_qa.tsv" \
    --tab_table \
    -q
else
  # Universal markers (faster)
  checkm tree \
    "$BINS_DIR" \
    "$OUTDIR/checkm_tree" \
    -x fa \
    -t "$THREADS" \
    2>&1 | tee "$OUTDIR/checkm.log"
  
  checkm tree_qa \
    "$OUTDIR/checkm_tree" \
    -o 2 \
    -f "$OUTDIR/checkm_qa.tsv" \
    --tab_table
fi

# Generate summary table with quality classification
python3 -c "
import sys
hq_total = mq_total = lq_total = bad_total = 0
with open('$OUTDIR/checkm_qa.tsv') as f:
    header = f.readline()
    print('bin\tcompleteness\tcontamination\tquality\thq_score\tmag_size')
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 5: continue
        bin_id = parts[0]
        try:
            comp = float(parts[1])
            cont = float(parts[2])
            size = int(parts[3]) if parts[3] else 0
        except:
            continue
        
        # MIMAG quality standards
        if comp >= 90 and cont < 5:
            quality = 'HQ'
        elif comp >= 50 and cont < 10:
            quality = 'MQ'
        elif comp >= 50 and cont < 5:
            quality = 'LQ'
        else:
            quality = 'Poor'
        
        # HQ score: completeness - 5*contamination
        hq_score = comp - 5 * cont
        
        if quality == 'HQ': hq_total += 1
        elif quality == 'MQ': mq_total += 1
        elif quality == 'LQ': lq_total += 1
        else: bad_total += 1
        
        print(f'{bin_id}\t{comp:.1f}\t{cont:.1f}\t{quality}\t{hq_score:.1f}\t{size}')

print(f'\\nSummary: HQ={hq_total}, MQ={mq_total}, LQ={lq_total}, Poor={bad_total}, Total={hq_total+mq_total+lq_total+bad_total}')
" > "$OUTDIR/checkm_summary.tsv"

cat "$OUTDIR/checkm_summary.tsv"

echo "─── CheckM Summary ───────────────────"
echo "  Total bins: $N_BINS"
echo "  Output: $OUTDIR/checkm_summary.tsv"
echo "═══════════════════════════════════════"
echo "✓ CheckM complete"
