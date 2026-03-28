# MAG 分析流程

**Metagenome-Assembled Genome (MAG) 组装与分箱分析流程 — Omics Case Study Week 第4天。**

## 🔬 流程概览

```
原始测序数据 (raw reads)
        │
        ▼
   fastp 质控过滤
        │
        ▼
   MEGAHIT 组装
        │
        ▼
   分箱 (Binning)
   MetaBAT2 / MaxBin2 / CONCOCT
        │
        ▼
   DAS Tool 去冗余整合
        │
        ▼
   CheckM 质量评估  ──────►  GTDB-Tk 物种注释
        │                              │
        ▼                              ▼
   Prokka 基因注释            CoverM 丰度估算
```

## 📋 工具链

| 步骤 | 工具 | 说明 |
|------|------|------|
| 1. 质控 | fastp | 过滤低质量 reads，去除接头 |
| 2. 组装 | MEGAHIT | 将过滤后 reads 组装为 contigs |
| 3a. 分箱 | MetaBAT2 | 基于覆盖度的分箱 |
| 3b. 分箱 | MaxBin2 | 基于标记基因的分箱 |
| 3c. 分箱 | CONCOCT | 基于覆盖度和组成的 NMF 分箱 |
| 4. 去冗余 | DAS Tool | 合并并去冗余多个分箱结果 |
| 5. 质量评估 | CheckM | 评估 Completeness 和 Contamination（MIMAG 标准） |
| 6. 物种注释 | GTDB-Tk | 基于基因组的分类学注释 |
| 7. 基因注释 | Prokka | 基因预测与功能注释 |
| 8. 丰度估算 | CoverM | MAG 在各样本中的相对丰度 |

## 📦 MIMAG 质量标准

| 质量等级 | Completeness（完整度） | Contamination（污染度） |
|----------|----------------------|------------------------|
| **High-quality (HQ)** | ≥ 90% | < 5% |
| **Medium-quality (MQ)** | ≥ 50% | < 10% |
| **Low-quality (LQ)** | ≥ 50% | < 5% |

> **注意**：High-quality MAG 还需要满足 23S、16S、5S rRNA 基因均存在且 tRNA 包含至少 18 种氨基酸。

## 🗂 目录结构

```
metagenomics-mag-pipeline/
├── README.md
├── Snakefile
├── .gitignore
├── config/
│   └── config.yaml              # 配置文件
├── bin/
│   ├── 01_qc.sh                 # 质控
│   ├── 02_assemble.sh           # MEGAHIT 组装
│   ├── 03a_metabat2.sh         # MetaBAT2 分箱
│   ├── 03b_maxbin2.sh          # MaxBin2 分箱
│   ├── 03c_concoct.sh          # CONCOCT 分箱
│   ├── 04_dastool.sh           # DAS Tool 去冗余
│   ├── 05_checkm.sh            # CheckM 质量评估
│   ├── 06_gtdbtk.sh            # GTDB-Tk 物种注释
│   ├── 07_prokka.sh            # Prokka 基因注释
│   └── 08_coverm.sh            # CoverM 丰度估算
├── scripts/
│   ├── summary_table.py        # 生成汇总表
│   └── plot_mag_quality.py     # 绘图脚本
├── envs/
│   └── environment.yaml        # Conda 环境配置
├── data/                        # 原始 reads（可用 symlink）
├── results/
│   ├── 01_qc/
│   ├── 02_assembly/
│   ├── 03_binning/
│   ├── 04_dastool/
│   ├── 05_checkm/
│   ├── 06_gtdbtk/
│   ├── 07_prokka/
│   └── 08_coverm/
└── logs/                        # Slurm/集群日志
```

## 🚀 安装

### 1. 克隆仓库

```bash
git clone https://github.com/nullvoid42/metagenomics-mag-pipeline.git
cd metagenomics-mag-pipeline
```

### 2. 创建 Conda 环境

```bash
conda env create -f envs/environment.yaml
conda activate mag-pipeline
```

### 3. 下载 GTDB-Tk 参考数据（约 5GB，需网络）

```bash
# 方法一：使用脚本（如果提供）
./bin/download_gtdbtk_refdata.sh

# 方法二：直接用 gtdbtk 命令
gtdbtk download_refdata --db_dir /path/to/gtdb_ref
```

> **HPC 服务器上安装示例**（连接到 `ssh wlab-weibin`）：
> ```bash
> git clone https://github.com/nullvoid42/metagenomics-mag-pipeline.git
> cd metagenomics-mag-pipeline
> conda env create -f envs/environment.yaml
> conda activate mag-pipeline
> gtdbtk download_refdata --db_dir /path/to/gtdb_ref
> ```

## ⚙️ 配置

编辑 `config/config.yaml`：

```yaml
project:
  name: "mag_pipeline"
  samples: ["sample_A", "sample_B"]   # 修改为你的样本名

data:
  reads_dir: "data/reads"
  adapters: ""                         # 自动检测，或提供 FASTA 文件

assembly:
  tool: "megahit"
  min_contig_len: 1500                 # 分箱最小 contig 长度

binners:
  metabat2: true
  maxbin2: true
  concoct: true

dastool:
  score_thresh: 0.5
  search_engine: "blast"

checkm:
  lineage_wf: true                     # 使用 lineage-specific marker set

gtdbtk:
  genome_path: "results/04_dastool/dereplicated_bins"
  out_path: "results/06_gtdbtk"

prokka:
  kingdom: "Bacteria"
  add_genes: true

coverm:
  method: "coverage"
  min_read_identity: 0.95
```

## 📥 测试数据下载

```bash
# 选项 A：CAMII SIM001 模拟土壤宏基因组（推荐，有 ground truth）
./bin/download_test_data.sh cami

# 选项 B：SRA 人类肠道宏基因组 ERR011100
./bin/download_test_data.sh sra

# 选项 C：使用你自己的数据 — 在 config/config.yaml 中设置路径
```

## ▶️ 运行

```bash
# 验证流程（dry run）
snakemake -n

# 本地运行（测试用，少量数据）
snakemake --use-conda -j 8

# HPC Slurm 集群提交
snakemake --use-conda \
  --cluster "sbatch -n 1 -t {params.runtime} --cpus-per-task {threads}" \
  -j 64

# 只运行特定步骤
snakemake --use-conda 05_checkm
```

## 📊 结果解读

### CheckM 质量报告

文件：`results/05_checkm/checkm_summary.tsv`

| bin | completeness | contamination | quality |
|-----|-------------|---------------|---------|
| bin.1 | 94.5 | 2.1 | HQ |
| bin.2 | 67.3 | 8.4 | MQ |

**如何看**：
- Completeness 越高越好（表示基因组越完整）
- Contamination 越低越好（表示混入的其他基因组越少）
- 优先使用 HQ MAG（≥90% complete, <5% contamination）进行下游分析
- MQ 和 LQ MAG 可用于功能基因分析，但不适合构建系统发育树

### GTDB-Tk 物种注释

文件：`results/06_gtdbtk/gtdbtk_summary.tsv`

| user_genome | classification | confidence |
|-------------|----------------|------------|
| bin.1 | d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria... | 0.95 |

**如何看**：
- `d__Bacteria` = Domain（域）
- `p__Proteobacteria` = Phylum（门）
- `c__Gammaproteobacteria` = Class（纲）
- `o__...;f__...;g__...;s__...` = Order/Family/Genus/Species
- confidence > 0.7 通常表示注释可信

### 汇总表

运行 `scripts/summary_table.py` 生成：
`results/summary_table.tsv`

包含：bin 名称、CheckM 质量、GTDB-Tk taxonomy、功能注释（Prokka）、丰度信息（CoverM）

### 图表

运行 `scripts/plot_mag_quality.py` 生成：

- `results/figures/mag_quality_barplot.pdf` — 各 MAG 的 Completeness vs Contamination
- `results/figures/taxonomy_piechart.pdf` — 门（Phylum）水平分类分布
- `results/figures/abundance_heatmap.pdf` — 样本 × MAG 丰度热图

## 📁 预期输出文件

```
results/
├── 01_qc/
│   └── *.clean.fastq.gz              # 质控后的 reads
├── 02_assembly/
│   └── megahit/
│       └── final.contigs.fa         # 组装结果
├── 03_binning/
│   ├── metabat2/
│   ├── maxbin2/
│   └── concoct/
├── 04_dastool/
│   └── dereplicated_bins/            # 去冗余后的 MAG
├── 05_checkm/
│   └── checkm_summary.tsv            # 质量评估结果
├── 06_gtdbtk/
│   └── gtdbtk_summary.tsv             # GTDB-Tk 分类学注释
├── 07_prokka/
│   └── *.faa                         # 蛋白序列
├── 08_coverm/
│   └── coverm_summary.tsv            # 丰度矩阵
└── figures/                          # 图表输出
```

## ❓ 常见问题

**Q: 组装后 contig 数量很多，但分箱效果差？**
A: 检查测序深度是否足够（建议 ≥ 1 Gbp per sample）。同时确认 `--min-contig-len` 参数是否合适（默认 1500）。

**Q: CheckM 报 `lineage-specific marker set not found`？**
A: 需要使用 `checkm lineage_wf` 或提前运行 `checkm test` 下载测试数据。也可手动指定 `--reduced-tree` 降低内存占用。

**Q: GTDB-Tk 运行很慢？**
A: GTDB-Tk 首次运行会下载 ~5GB 参考数据。后续可设置 `GTDBTK_DATA_PATH` 环境变量指向缓存目录加速。也可在 `config.yaml` 中加 `--force` 重新运行。

**Q: DAS Tool 报 `No bins above score threshold`？**
A: 检查各分箱工具（MetaBAT2/MaxBin2/CONCOCT）是否成功生成了 bins。确认 score threshold 不要设得太高（默认 0.5）。

**Q: Prokka 注释结果太少？**
A: 对于高污染或低质量 MAG，Prokka 可能检测不到足够基因。先过滤掉 LQ MAG，再单独运行 Prokka。

**Q: 如何只分析部分样本？**
A: 在 `config.yaml` 的 `samples` 列表中只保留需要的样本名，重新运行即可。

## 📚 引用

- Sealy et al. (2024). GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. *Bioinformatics*.
- Kang et al. (2019). MetaBAT2: an efficient binning algorithm for large metagenome assemblies. *BMC Bioinformatics*.
- Wu et al. (2016). MaxBin 2.0: an automated binning algorithm for recovering genomes from microbial communities. *Microbiome*.
- Alneberg et al. (2014). CONCOCT: clustering contigs on coverage and composition. *PLoS ONE*.
- Sieber et al. (2018). Recovery of genomes from metagenomes via a dereplication strategy. *Nature Microbiology*.
- Parks et al. (2020). GTDB: an updated census of the bacterial Tree of Life. *F1000Research*.
- Bowers et al. (2017). Minimum information about a metagenome-assembled genome (MIMAG). *Nature Biotechnology*.

## 📝 License

MIT — freely usable, cite if helpful.
