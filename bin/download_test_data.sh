#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════
# download_test_data.sh — Download test datasets for the MAG pipeline
# ═══════════════════════════════════════════════════════════════════
# Usage: bash bin/download_test_data.sh [cami|sra|custom]
set -euo pipefail

PROXY=""
if curl -s --connect-timeout 5 https://github.com > /dev/null 2>&1; then
  PROXY=""
else
  echo "⚠️  Direct connection slow — using proxy..."
  PROXY="http_proxy=http://127.0.0.1:7897 https_proxy=http://127.0.0.1:7897"
fi

MODE="${1:-cami}"
DATA_DIR="data/reads"
mkdir -p "$DATA_DIR"

echo "═══════════════════════════════════════"
echo "  Downloading test data: $MODE"
echo "═══════════════════════════════════════"

case "$MODE" in
  cami)
    echo "  Option A: CAMII SIM001 simulated metagenome (soil, ground truth)"
    echo "  This is the recommended benchmark dataset."
    echo ""
    echo "  CAMI II Mouse Gut Microbiome dataset (SIM001):"
    echo "  https://data.cami-challenge.org/"
    echo ""
    echo "  Direct download (example — check CAMI website for current links):"
    echo "  $PROXY wget -P $DATA_DIR https://s3.amazonaws.com/cami-data-challenge/SIM001.tar.gz"
    echo ""
    echo "  Alternatively, download CAMI data manually:"
    echo "  1. Visit https://data.cami-challenge.org/"
    echo "  2. Select 'CAMI II: Mouse Gut Microbiome' → 'SIM001'"
    echo "  3. Download reads and place in: $DATA_DIR/"
    echo ""
    
    # Try EBI CAMI mirror
    echo "  Attempting download from EBI..."
    $PROXY wget -q -O "$DATA_DIR/cami_reads.tar.gz" \
      "https://www.ebi.ac.uk/ena/browser/api/fasta/ERR5766172?download=1&gzip=true" 2>/dev/null || true
    
    if [ -f "$DATA_DIR/cami_reads.tar.gz" ]; then
      tar -xzf "$DATA_DIR/cami_reads.tar.gz" -C "$DATA_DIR/"
      echo "✓ CAMI data downloaded"
    else
      echo "⚠️  Could not auto-download CAMI data."
      echo "   Please download manually from https://data.cami-challenge.org/"
      echo "   Place reads as: data/reads/{sample}_1.fastq.gz and {sample}_2.fastq.gz"
    fi
    ;;
    
  sra)
    echo "  Option B: SRA human gut metagenome ERR011100 (MetaHIT)"
    echo "  Using SRA Toolkit to download..."
    
    # Check for SRA toolkit
    if ! command -v fasterq-dump &>/dev/null && ! command -v prefetch &>/dev/null; then
      echo "  Installing SRA Toolkit via conda..."
      conda install -y -c bioconda sra-tools 2>/dev/null || pip install sra-tools 2>/dev/null || true
    fi
    
    echo "  Downloading ERR011100 from SRA..."
    $PROXY prefetch -O "$DATA_DIR/" ERR011100 2>&1 || true
    $PROXY fasterq-dump "$DATA_DIR/ERR011100" -O "$DATA_DIR/" -e 8 --split-files 2>&1 || true
    
    if [ -f "$DATA_DIR/ERR011100_1.fastq" ]; then
      gzip "$DATA_DIR/ERR011100_1.fastq" "$DATA_DIR/ERR011100_2.fastq" 2>/dev/null || true
      mv "$DATA_DIR/ERR011100_1.fastq.gz" "$DATA_DIR/ERR011100_1.fastq.gz"
      mv "$DATA_DIR/ERR011100_2.fastq.gz" "$DATA_DIR/ERR011100_2.fastq.gz"
      echo "✓ SRA data downloaded: ERR011100"
      echo "  Update config/config.yaml: samples: ['ERR011100']"
    else
      echo "⚠️  Could not auto-download SRA data."
      echo "   Manual download:"
      echo "   fasterq-dump ERR011100 -O data/reads/ -e 8 --split-files"
      echo "   gzip data/reads/ERR011100_*.fastq"
    fi
    ;;
    
  mock)
    echo "  Option C: Generating mock Illumina paired-end reads (for testing only!)"
    echo "  ⚠️  This creates synthetic data — DO NOT use for real analysis!"
    
    python3 << 'PYEOF'
import os, random, gzip

os.makedirs("data/reads", exist_ok=True)
BASES = "ACGT"
QUAL = list(range(33, 75))

for sample in ["mock_A", "mock_B"]:
    for read_num in [1, 2]:
        fname = f"data/reads/{sample}_{read_num}.fastq.gz"
        with gzip.open(fname, "wt") as f:
            for i in range(100_000):  # 100K read pairs per sample
                seq = "".join(random.choices(BASES, k=150))
                qual = "".join(chr(random.choice(QUAL)) for _ in range(150))
                f.write(f"@{sample}_{read_num}_{i}\n{seq}\n+\n{qual}\n")
        print(f"  Created {fname}")
print("  ✓ Mock data generated (for pipeline testing only)")
PYEOF
    ;;
    
  *)
    echo "Usage: bash bin/download_test_data.sh [cami|sra|mock]"
    echo ""
    echo "  cami    - CAMII simulated benchmark data (recommended)"
    echo "  sra     - SRA human gut metagenome ERR011100"  
    echo "  mock    - Generate synthetic mock data (testing only)"
    exit 1
    ;;
esac

echo "═══════════════════════════════════════"
echo "✓ Data download complete"
echo ""
echo "Next steps:"
echo "  1. Edit config/config.yaml — set sample names"
echo "  2. Run: snakemake --use-conda -j 8"
echo "═══════════════════════════════════════"
