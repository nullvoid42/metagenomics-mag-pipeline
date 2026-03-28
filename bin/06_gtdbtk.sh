#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════
# 06_gtdbtk.sh — GTDB-Tk taxonomy classification
# ═══════════════════════════════════════════════════════════════════
# Usage: bash bin/06_gtdbtk.sh <bins_dir> <outdir> <threads>
set -euo pipefail

BINS_DIR="$1"; OUTDIR="$2"; THREADS="${3:-32}"

mkdir -p "$OUTDIR"

echo "═══════════════════════════════════════"
echo "  GTDB-Tk taxonomy classification"
echo "  Bins dir: $BINS_DIR"
echo "═══════════════════════════════════════"

# Check GTDBTK_DATA_PATH
if [ -z "${GTDBTK_DATA_PATH:-}" ]; then
  echo "WARNING: GTDBTK_DATA_PATH not set. Trying default locations..."
  if [ -d "/db2/gtdbtk" ]; then
    GTDBTK_DATA_PATH="/db2/gtdbtk"
    export GTDBTK_DATA_PATH
    echo "  Found GTDB-Tk data at: $GTDBTK_DATA_PATH"
  else
    echo "ERROR: GTDB-Tk reference data not found."
    echo "  Download with: gtdbtk download_refdata --db_dir <path>"
    echo "  Or set: export GTDBTK_DATA_PATH=/path/to/gtdb_ref"
    exit 1
  fi
fi

N_BINS=$(find "$BINS_DIR" -name "*.fa" -o -name "*.fasta" 2>/dev/null | wc -l | tr -d ' ')
echo "  Found $N_BINS bins to classify"

# ── gtdbtk classify_wf ────────────────────────────────────────────
gtdbtk classify_wf \
  --genome_dir "$BINS_DIR" \
  --out_dir "$OUTDIR" \
  --extension fa \
  --cpus "$THREADS" \
  --force \
  2>&1 | tee "$OUTDIR/gtdbtk.log"

# ── Summarize results ────────────────────────────────────────────
SUMMARY_TSV="$OUTDIR/gtdbtk_summary.tsv"

# Extract from gtdbtk summary
if [ -f "$OUTDIR/gtdbtk.summary.tsv" ]; then
  tail -n +2 "$OUTDIR/gtdbtk.summary.tsv" | \
    awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > "$SUMMARY_TSV"
fi

# Add taxonomy column
python3 -c "
import sys
print('user_genome\tclassification\tnote\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies')
with open('$OUTDIR/gtdbtk.bac120.summary.tsv' if '$OUTDIR/gtdbtk.bac120.summary.tsv' else '$OUTDIR/gtdbtk.summary.tsv', 'r') as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 3: continue
        genome = parts[0]
        classification = parts[1] if len(parts) > 1 else 'Unknown'
        note = parts[2] if len(parts) > 2 else ''
        
        # Parse taxonomy
        taxonomy = {}
        for rank in ['d__','p__','c__','o__','f__','g__','s__']:
            if rank in classification:
                val = classification.split(rank)[1].split(';')[0] if ';' in classification else ''
                taxonomy[rank] = val
            else:
                taxonomy[rank] = ''
        
        d = taxonomy.get('d__','')
        p = taxonomy.get('p__','')
        c = taxonomy.get('c__','')
        o = taxonomy.get('o__','')
        f = taxonomy.get('f__','')
        g = taxonomy.get('g__','')
        s = taxonomy.get('s__','')
        
        print(f'{genome}\t{classification}\t{note}\t{d}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}')
" > "$SUMMARY_TSV"

echo "─── GTDB-Tk Summary ───────────────────"
echo "  Classified: $N_BINS bins"
echo "  Output: $SUMMARY_TSV"
echo "═══════════════════════════════════════"
echo "✓ GTDB-Tk complete"
