# Snakemake workflow for MAG (Metagenome-Assembled Genome) analysis
# Day 4: Omics Case Study Week

import os
import glob

# ─── Configuration ───────────────────────────────────────────────────────────
configfile: "config/config.yaml"

PROJECT_NAME = config["project"]["name"]
SAMPLES = config["project"]["samples"]
READS_DIR = config["data"]["reads_dir"]
ADAPTERS = config["data"].get("adapters", "")
THREADS = config["project"].get("threads", 16)

# Derived paths
QC_DIR = "results/01_qc"
ASSEMBLY_DIR = "results/02_assembly"
BINNING_DIR = "results/03_binning"
DASTOOL_DIR = "results/04_dastool"
CHECKM_DIR = "results/05_checkm"
GTDBTK_DIR = "results/06_gtdbtk"
PROKKA_DIR = "results/07_prokka"
COVERM_DIR = "results/08_coverm"
SCRIPTS_DIR = "scripts"
FIGURES_DIR = "results/figures"

# Quality thresholds (MIMAG)
HQ_COMPLETENESS = 90
HQ_CONTAMINATION = 5
MQ_COMPLETENESS = 50
MQ_CONTAMINATION = 10

# ─── Helper functions ────────────────────────────────────────────────────────
def get_sample_reads(wildcards):
    """Return list of read files for a sample."""
    sample = wildcards.sample
    r1 = os.path.join(READS_DIR, f"{sample}_1.fastq.gz")
    r2 = os.path.join(READS_DIR, f"{sample}_2.fastq.gz")
    return {"r1": r1, "r2": r2}

def get_dastool_bins():
    """Return list of dereplicated bin files."""
    bins = glob.glob(os.path.join(DASTOOL_DIR, "dereplicated_bins", "*.fa"))
    return bins

def get_all_binners():
    """Return list of all binner output directories."""
    binners = []
    if config["binners"].get("metabat2", True):
        binners.append(os.path.join(BINNING_DIR, "metabat2"))
    if config["binners"].get("maxbin2", True):
        binners.append(os.path.join(BINNING_DIR, "maxbin2"))
    if config["binners"].get("concoct", True):
        binners.append(os.path.join(BINNING_DIR, "concoct"))
    return binners

# ─── Download reference data (GTDB-Tk) ──────────────────────────────────────
rule download_gtdbtk_ref:
    output:
        touch(os.path.join(ASSEMBLY_DIR, ".gtdbtk_ref_done"))
    shell:
        """
        if [ -z "$GTDBTK_DATA_PATH" ]; then
            echo "ERROR: Set GTDBTK_DATA_PATH environment variable first."
            echo "Download: gtdbtk download_refdata --db_dir $GTDBTK_DATA_PATH"
            exit 1
        fi
        echo "GTDB-Tk reference data should already be at: $GTDBTK_DATA_PATH"
        """

# ─── Step 1: Quality control with fastp ────────────────────────────────────
rule qc_reads:
    input:
        unpack(get_sample_reads)
    output:
        r1 = os.path.join(QC_DIR, "{sample}_1.filtered.fastq.gz"),
        r2 = os.path.join(QC_DIR, "{sample}_2.filtered.fastq.gz"),
        json = os.path.join(QC_DIR, "{sample}.fastp.json"),
        html = os.path.join(QC_DIR, "{sample}.fastp.html")
    log:
        os.path.join("logs", "01_qc", "{sample}.log")
    params:
        adapters = ADAPTERS,
        outdir = QC_DIR
    shell:
        """
        mkdir -p {params.outdir} {params.outdir}/unpaired
        bash bin/01_qc.sh {input.r1} {input.r2} {params.outdir} {wildcards.sample} "{params.adapters}" 2>&1 | tee {log}
        """

# ─── Step 2: Assembly with MEGAHIT ─────────────────────────────────────────
rule assemble_megahit:
    input:
        r1 = os.path.join(QC_DIR, "{sample}_1.filtered.fastq.gz"),
        r2 = os.path.join(QC_DIR, "{sample}_2.filtered.fastq.gz")
    output:
        contigs = os.path.join(ASSEMBLY_DIR, "{sample}", "{sample}.contigs.fa"),
        done = os.path.join(ASSEMBLY_DIR, "{sample}", "done")
    log:
        os.path.join("logs", "02_assembly", "{sample}.log")
    params:
        outdir = os.path.join(ASSEMBLY_DIR, "{sample}"),
        min_len = config["assembly"].get("min_contig_len", 1500)
    shell:
        """
        mkdir -p {params.outdir} logs/02_assembly
        bash bin/02_assemble.sh {input.r1} {input.r2} {params.outdir} {wildcards.sample} {params.min_len} 2>&1 | tee {log}
        touch {output.done}
        """

# ─── Step 3a: MetaBAT2 binning ──────────────────────────────────────────────
rule binning_metabat2:
    input:
        contigs = os.path.join(ASSEMBLY_DIR, "{sample}", "{sample}.contigs.fa"),
        done = os.path.join(ASSEMBLY_DIR, "{sample}", "done")
    output:
        touch(os.path.join(BINNING_DIR, "metabat2", "{sample}", ".done"))
    log:
        os.path.join("logs", "03a_metabat2", "{sample}.log")
    params:
        outdir = os.path.join(BINNING_DIR, "metabat2", "{sample}"),
        min_len = config["assembly"].get("min_contig_len", 1500)
    shell:
        """
        mkdir -p {params.outdir} logs/03a_metabat2 logs/03b_maxbin2 logs/03c_concoct
        bash bin/03a_metabat2.sh {input.contigs} {params.outdir} {params.min_len} 2>&1 | tee {log}
        """

# ─── Step 3b: MaxBin2 binning ───────────────────────────────────────────────
rule binning_maxbin2:
    input:
        contigs = os.path.join(ASSEMBLY_DIR, "{sample}", "{sample}.contigs.fa"),
        done = os.path.join(ASSEMBLY_DIR, "{sample}", "done")
    output:
        touch(os.path.join(BINNING_DIR, "maxbin2", "{sample}", ".done"))
    log:
        os.path.join("logs", "03b_maxbin2", "{sample}.log")
    params:
        outdir = os.path.join(BINNING_DIR, "maxbin2", "{sample}"),
        min_len = config["assembly"].get("min_contig_len", 1500)
    shell:
        """
        mkdir -p {params.outdir}
        bash bin/03b_maxbin2.sh {input.contigs} {params.outdir} {params.min_len} 2>&1 | tee {log}
        """

# ─── Step 3c: CONCOCT binning ────────────────────────────────────────────────
rule binning_concoct:
    input:
        contigs = os.path.join(ASSEMBLY_DIR, "{sample}", "{sample}.contigs.fa"),
        done = os.path.join(ASSEMBLY_DIR, "{sample}", "done")
    output:
        touch(os.path.join(BINNING_DIR, "concoct", "{sample}", ".done"))
    log:
        os.path.join("logs", "03c_concoct", "{sample}.log")
    params:
        outdir = os.path.join(BINNING_DIR, "concoct", "{sample}"),
        min_len = config["assembly"].get("min_contig_len", 1500)
    shell:
        """
        mkdir -p {params.outdir}
        bash bin/03c_concoct.sh {input.contigs} {params.outdir} {params.min_len} 2>&1 | tee {log}
        """

# ─── Collect all binning results ────────────────────────────────────────────
rule collect_binning:
    input:
        expand(os.path.join(BINNING_DIR, "metabat2", "{sample}", ".done"), sample=SAMPLES),
        expand(os.path.join(BINNING_DIR, "maxbin2", "{sample}", ".done"), sample=SAMPLES),
        expand(os.path.join(BINNING_DIR, "concoct", "{sample}", ".done"), sample=SAMPLES)
    output:
        touch(os.path.join(BINNING_DIR, ".all_binners_done"))
    shell:
        "echo 'All binners completed'"

# ─── Step 4: DAS Tool dereplication ─────────────────────────────────────────
rule dastool_dereplicate:
    input:
        bin_done = os.path.join(BINNING_DIR, ".all_binners_done")
    output:
        directory(os.path.join(DASTOOL_DIR, "dereplicated_bins")),
        os.path.join(DASTOOL_DIR, "DASTool_summary.tsv")
    log:
        os.path.join("logs", "04_dastool", "dastool.log")
    params:
        metabat2_dir = os.path.join(BINNING_DIR, "metabat2"),
        maxbin2_dir = os.path.join(BINNING_DIR, "maxbin2"),
        concoct_dir = os.path.join(BINNING_DIR, "concoct"),
        outdir = DASTOOL_DIR,
        samples = SAMPLES,
        min_score = config["dastool"].get("score_thresh", 0.5)
    shell:
        """
        mkdir -p {params.outdir}/dereplicated_bins logs/04_dastool
        bash bin/04_dastool.sh {params.metabat2_dir} {params.maxbin2_dir} {params.concoct_dir} \\
            {params.outdir} {params.min_score} {params.samples} 2>&1 | tee {log}
        """

# ─── Step 5: CheckM quality assessment ───────────────────────────────────────
rule checkm_assess:
    input:
        bins_dir = os.path.join(DASTOOL_DIR, "dereplicated_bins")
    output:
        os.path.join(CHECKM_DIR, "checkm_summary.tsv")
    log:
        os.path.join("logs", "05_checkm", "checkm.log")
    params:
        outdir = CHECKM_DIR,
        lineage = config["checkm"].get("lineage_wf", True),
        threads = THREADS
    shell:
        """
        mkdir -p {params.outdir} logs/05_checkm
        bash bin/05_checkm.sh {input.bins_dir} {params.outdir} {params.lineage} {params.threads} 2>&1 | tee {log}
        """

# ─── Step 6: GTDB-Tk taxonomy classification ───────────────────────────────
rule gtdbtk_classify:
    input:
        bins_dir = os.path.join(DASTOOL_DIR, "dereplicated_bins")
    output:
        os.path.join(GTDBTK_DIR, "gtdbtk_summary.tsv")
    log:
        os.path.join("logs", "06_gtdbtk", "gtdbtk.log")
    params:
        outdir = GTDBTK_DIR,
        threads = THREADS
    shell:
        """
        mkdir -p {params.outdir} logs/06_gtdbtk
        bash bin/06_gtdbtk.sh {input.bins_dir} {params.outdir} {params.threads} 2>&1 | tee {log}
        """

# ─── Step 7: Prokka gene annotation ─────────────────────────────────────────
rule prokka_annotate:
    input:
        os.path.join(DASTOOL_DIR, "dereplicated_bins")
    output:
        directory(os.path.join(PROKKA_DIR, "annotations"))
    log:
        os.path.join("logs", "07_prokka", "prokka.log")
    params:
        outdir = os.path.join(PROKKA_DIR, "annotations"),
        kingdom = config["prokka"].get("kingdom", "Bacteria"),
        threads = THREADS
    shell:
        """
        mkdir -p {params.outdir} logs/07_prokka
        bash bin/07_prokka.sh {input} {params.outdir} {params.kingdom} {params.threads} 2>&1 | tee {log}
        """

# ─── Step 8: CoverM abundance estimation ───────────────────────────────────
rule coverm_abundance:
    input:
        bins = os.path.join(DASTOOL_DIR, "dereplicated_bins"),
        r1 = expand(os.path.join(QC_DIR, "{sample}_1.filtered.fastq.gz"), sample=SAMPLES),
        r2 = expand(os.path.join(QC_DIR, "{sample}_2.filtered.fastq.gz"), sample=SAMPLES)
    output:
        os.path.join(COVERM_DIR, "coverm_summary.tsv")
    log:
        os.path.join("logs", "08_coverm", "coverm.log")
    params:
        outdir = COVERM_DIR,
        threads = THREADS,
        samples = SAMPLES
    shell:
        """
        mkdir -p {params.outdir} logs/08_coverm
        bash bin/08_coverm.sh {input.bins} {params.outdir} {params.threads} {params.samples} 2>&1 | tee {log}
        """

# ─── Generate summary table ──────────────────────────────────────────────────
rule summary_table:
    input:
        checkm = os.path.join(CHECKM_DIR, "checkm_summary.tsv"),
        gtdbtk = os.path.join(GTDBTK_DIR, "gtdbtk_summary.tsv"),
        coverm = os.path.join(COVERM_DIR, "coverm_summary.tsv")
    output:
        os.path.join("results", "summary_table.tsv")
    shell:
        """
        python3 scripts/summary_table.py \\
            --checkm {input.checkm} \\
            --gtdbtk {input.gtdbtk} \\
            --coverm {input.coverm} \\
            --output {output}
        """

# ─── Generate publication figures ───────────────────────────────────────────
rule plot_figures:
    input:
        os.path.join("results", "summary_table.tsv")
    output:
        directory(FIGURES_DIR)
    shell:
        """
        mkdir -p {output}
        python3 scripts/plot_mag_quality.py \\
            --input {input} \\
            --outdir {output}
        """

# ─── Full pipeline ────────────────────────────────────────────────────────────
rule all:
    default_target: True
    input:
        os.path.join("results", "summary_table.tsv"),
        FIGURES_DIR

# ─── Cluster configuration ────────────────────────────────────────────────────
# Use with: snakemake --cluster "sbatch -n 1 -t {params.runtime} --cpus-per-task {threads}"
# Example cluster config passed via --cluster-config
cluster_config = {
    "qc_reads": {"runtime": "02:00:00", "threads": 8},
    "assemble_megahit": {"runtime": "12:00:00", "threads": 32},
    "binning_metabat2": {"runtime": "08:00:00", "threads": 16},
    "binning_maxbin2": {"runtime": "08:00:00", "threads": 16},
    "binning_concoct": {"runtime": "08:00:00", "threads": 16},
    "dastool_dereplicate": {"runtime": "04:00:00", "threads": 16},
    "checkm_assess": {"runtime": "08:00:00", "threads": 16},
    "gtdbtk_classify": {"runtime": "12:00:00", "threads": 32},
    "prokka_annotate": {"runtime": "06:00:00", "threads": 16},
    "coverm_abundance": {"runtime": "04:00:00", "threads": 16},
    "summary_table": {"runtime": "01:00:00", "threads": 4},
    "plot_figures": {"runtime": "01:00:00", "threads": 4},
}

# ─── Report for MultiQC ───────────────────────────────────────────────────────
ruleorder: qc_reads > assemble_megahit

onsuccess:
    print("\n✅ Pipeline completed successfully!")
    print("📊 Results: results/summary_table.tsv")
    print("📈 Figures: results/figures/")
