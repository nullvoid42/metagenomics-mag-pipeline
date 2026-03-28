#!/usr/bin/env python3
"""
summary_table.py — Aggregate CheckM, GTDB-Tk, and CoverM results into a final TSV table.
Outputs: results/summary_table.tsv

Usage:
    python3 scripts/summary_table.py \\
        --checkm results/05_checkm/checkm_summary.tsv \\
        --gtdbtk results/06_gtdbtk/gtdbtk_summary.tsv \\
        --coverm results/08_coverm/coverm_summary.tsv \\
        --output results/summary_table.tsv
"""

import argparse
import os
import sys
from collections import defaultdict

def parse_checkm(path):
    """Parse CheckM summary: bin, completeness, contamination, quality, hq_score, mag_size"""
    bins = {}
    if not os.path.exists(path):
        print(f"WARNING: CheckM file not found: {path}", file=sys.stderr)
        return bins
    with open(path) as f:
        header = f.readline()
        for line in f:
            line = line.strip()
            if not line or line.startswith("Summary"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            bin_id = os.path.basename(parts[0]).replace(".fa", "").replace(".fasta", "")
            bins[bin_id] = {
                "completeness": float(parts[1]) if parts[1] else 0.0,
                "contamination": float(parts[2]) if parts[2] else 100.0,
                "quality": parts[3],
                "hq_score": float(parts[4]) if parts[4] else 0.0,
                "mag_size": int(parts[5]) if parts[5] else 0,
            }
    return bins

def parse_gtdbtk(path):
    """Parse GTDB-Tk summary: user_genome, classification, phylum, genus, species"""
    taxonomy = {}
    if not os.path.exists(path):
        print(f"WARNING: GTDB-Tk file not found: {path}", file=sys.stderr)
        return taxonomy
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            genome = os.path.basename(parts[0]).replace(".fa", "").replace(".fasta", "")
            taxonomy[genome] = {
                "classification": parts[1],
                "domain": parts[3] if len(parts) > 3 else "",
                "phylum": parts[4] if len(parts) > 4 else "",
                "class": parts[5] if len(parts) > 5 else "",
                "order": parts[6] if len(parts) > 6 else "",
                "family": parts[7] if len(parts) > 7 else "",
                "genus": parts[8] if len(parts) > 8 else "",
                "species": parts[9] if len(parts) > 9 else "",
            }
    return taxonomy

def parse_coverm(path):
    """Parse CoverM summary: mag, genus, abundance, coverage, length, count"""
    abundance = {}
    if not os.path.exists(path):
        print(f"WARNING: CoverM file not found: {path}", file=sys.stderr)
        return abundance
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            mag = os.path.basename(parts[0]).replace(".fa", "").replace(".fasta", "")
            abundance[mag] = {
                "genus": parts[1],
                "abundance": float(parts[2]) if parts[2] else 0.0,
                "coverage": float(parts[3]) if len(parts) > 3 and parts[3] else 0.0,
                "length": int(parts[4]) if len(parts) > 4 and parts[4] else 0,
                "count": int(parts[5]) if len(parts) > 5 and parts[5] else 0,
            }
    return abundance

def main():
    parser = argparse.ArgumentParser(description="Aggregate MAG analysis results")
    parser.add_argument("--checkm", required=True, help="Path to CheckM summary TSV")
    parser.add_argument("--gtdbtk", required=True, help="Path to GTDB-Tk summary TSV")
    parser.add_argument("--coverm", required=True, help="Path to CoverM summary TSV")
    parser.add_argument("--output", default="results/summary_table.tsv", help="Output TSV")
    args = parser.parse_args()

    checkm = parse_checkm(args.checkm)
    gtdbtk = parse_gtdbtk(args.gtdbtk)
    coverm = parse_coverm(args.coverm)

    # Collect all MAG IDs
    all_mags = sorted(set(checkm.keys()) | set(gtdbtk.keys()) | set(coverm.keys()))

    # Header
    cols = [
        "mag_id", "quality", "completeness", "contamination", "hq_score",
        "mag_size_bp", "coverage", "abundance",
        "taxonomy_full", "domain", "phylum", "class", "order", "family", "genus", "species"
    ]

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w") as out:
        out.write("\t".join(cols) + "\n")
        for mag in all_mags:
            c = checkm.get(mag, {})
            g = gtdbtk.get(mag, {})
            v = coverm.get(mag, {})

            row = {
                "mag_id": mag,
                "quality": c.get("quality", "Unknown"),
                "completeness": c.get("completeness", ""),
                "contamination": c.get("contamination", ""),
                "hq_score": c.get("hq_score", ""),
                "mag_size_bp": c.get("mag_size", v.get("length", "")),
                "coverage": v.get("coverage", ""),
                "abundance": v.get("abundance", ""),
                "taxonomy_full": g.get("classification", "Unknown"),
                "domain": g.get("domain", ""),
                "phylum": g.get("phylum", ""),
                "class": g.get("class", ""),
                "order": g.get("order", ""),
                "family": g.get("family", ""),
                "genus": g.get("genus", ""),
                "species": g.get("species", ""),
            }
            out.write("\t".join(str(row[c]) for c in cols) + "\n")

    # Print summary stats
    hq = sum(1 for m in checkm.values() if m.get("quality") == "HQ")
    mq = sum(1 for m in checkm.values() if m.get("quality") == "MQ")
    lq = sum(1 for m in checkm.values() if m.get("quality") == "LQ")
    poor = sum(1 for m in checkm.values() if m.get("quality") == "Poor")
    total = len(checkm)

    print(f"\n{'='*50}")
    print(f"  MAG Summary Table: {args.output}")
    print(f"  Total MAGs: {total}")
    print(f"  High Quality (HQ, ≥90% complete, <5% contam): {hq}")
    print(f"  Medium Quality (MQ, ≥50% complete, <10% contam): {mq}")
    print(f"  Low Quality (LQ, ≥50% complete, <5% contam): {lq}")
    print(f"  Poor: {poor}")
    print(f"{'='*50}")

if __name__ == "__main__":
    main()
