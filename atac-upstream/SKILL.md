---
name: ATAC-seq Upstream Analysis
description: >
  This skill should be used when the user asks to "run ATAC-seq upstream analysis",
  "ATAC-seq preprocessing", "ATAC-seq QC", "ATAC-seq peak calling",
  "nf-core atacseq", "ENCODE ATAC-seq pipeline", "ATAC-seq比对",
  "ATAC-seq上游分析", "染色质开放性分析", "run ATAC-seq pipeline",
  "ATAC-seq alignment", "MACS2 peak calling", or mentions chromatin
  accessibility preprocessing workflows. Covers the full upstream pipeline from
  raw FASTQ through alignment, deduplication, peak calling, and QC reporting.
version: "1.0"
---

# ATAC-seq Upstream Analysis

This skill guides the complete upstream workflow for Assay for Transposase-Accessible Chromatin using sequencing (ATAC-seq). It covers every step from raw FASTQ files to final peak calls and QC reports, providing decision criteria for pipeline selection, environment setup, and quality gate thresholds.

## Purpose and Scope (目的与范围)

This skill covers the **upstream** portion of ATAC-seq analysis only. The pipeline stages included are:

1. **Raw data QC** (FastQC, adapter trimming)
2. **Read alignment** (Bowtie2 / BWA to reference genome)
3. **Duplicate removal** (Picard / SAMtools markdup)
4. **Mitochondrial read filtering**
5. **Peak calling** (MACS2, optionally with IDR reproducibility)
6. **Signal track generation** (bigWig coverage files)
7. **Comprehensive QC reporting** (MultiQC aggregation)

For downstream visualization, differential accessibility, motif analysis, and chromosome-level inspection, refer to the `atac-chromosome` skill.

## Environment Setup (环境准备)

Activate the pre-configured `upstream` mamba environment and verify all required tools:

```bash
conda activate upstream
```

Verify tool availability:

```bash
nextflow -version        # expect >= 25.10.4
caper --version          # expect >= 2.3.2
samtools --version
macs2 --version
multiqc --version
fastqc --version
bedtools --version
java -version             # expect openjdk 17
```

If any tool is missing or reports an unexpected version, reinstall it with `mamba install -n upstream <tool>`.

**Pipeline source code** is available locally at these paths:

| Pipeline | Local Path |
|----------|-----------|
| nf-core/atacseq | `/home/user/PanChongshi/Repo/nf-core-atacseq` |
| ENCODE-DCC atac-seq-pipeline | `/home/user/PanChongshi/Repo/encode-atac-seq-pipeline` |
| PEPATAC | `/home/user/PanChongshi/Repo/pepatac` |

## Pipeline Selection Decision Tree (流程选择指南)

Choose a pipeline based on the table below. The decision hinges on compliance requirements, compute environment, and project scale.

| Criterion | nf-core/atacseq | ENCODE-DCC | PEPATAC |
|-----------|----------------|------------|---------|
| **Best for** | General-purpose ATAC-seq | ENCODE-compliant production | Flexible research projects |
| **Workflow language** | Nextflow (DSL2) | WDL | Python / Looper |
| **Orchestration** | Nextflow executor | Caper / Cromwell | Looper (serial or cluster) |
| **Container support** | Docker, Conda, Singularity | Docker | Docker, Conda |
| **QC depth** | Good (MultiQC + built-in) | Excellent (ENCODE gold-standard) | Moderate |
| **IDR reproducibility** | Yes (optional) | Yes (built-in) | No |
| **Multi-sample aggregation** | Yes | Yes | Limited |
| **Ease of customization** | High (modular DSL2) | Medium (JSON configs) | High (PEP format) |
| **Regulatory compliance** | None specific | ENCODE DAC standards | None specific |
| **Learning curve** | Medium | Steep | Low |

**Quick decision rules:**

- Publishing to ENCODE or need ENCODE-compliant QC metrics? Use **ENCODE-DCC**.
- Running on an HPC cluster with SLURM/Singularity and want modular reproducibility? Use **nf-core/atacseq**.
- Prototyping with few samples, want minimal config overhead, or prefer Python ecosystem? Use **PEPATAC**.

## Quick Start: nf-core/atacseq

This is the recommended default for most projects. nf-core/atacseq provides a modular, containerized Nextflow workflow with strong community support.

**1. Prepare the samplesheet (CSV format):**

```csv
sample,fastq_1,fastq_2
S1,/data/S1_R1.fastq.gz,/data/S1_R2.fastq.gz
S2,/data/S2_R1.fastq.gz,/data/S2_R2.fastq.gz
```

**2. Build the genome index (first time only):**

```bash
nextflow run nf-core/atacseq \
  -profile docker \
  --genome GRCh38 \
  --build_index \
  --outdir /data/genome_index
```

**3. Run the full pipeline:**

```bash
nextflow run nf-core/atacseq \
  --input samplesheet.csv \
  --genome GRCh38 \
  --skip_ataqv false \
  --outdir results_atac \
  -profile docker
```

Replace `-profile docker` with `-profile singularity` on HPC clusters. Add `-resume` to recover from interrupted runs.

## Quick Start: ENCODE-DCC atac-seq-pipeline

Use this pipeline when ENCODE compliance is required or when gold-standard QC metrics are critical.

**1. Create the input JSON file:**

```json
{
  "atac.pipeline_type": "tf",
  "atac.genome_tsv": "https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v2/hg38.tsv",
  "atac.fastqs_R1": ["/data/S1_R1.fastq.gz"],
  "atac.fastqs_R2": ["/data/S1_R2.fastq.gz"],
  "atac.adapters": ["CTGTCTCTTATACACATCT"],
  "atac.qc_report.json": "qc_report.json"
}
```

**2. Run with caper:**

```bash
caper run atac.wdl \
  -i input.json \
  -b google -c caper.conf
```

Replace `-b google` with `-b local` for local execution or `-b slurm` for HPC clusters. The WDL file (`atac.wdl`) is located in the local clone at `/home/user/PanChongshi/Repo/encode-atac-seq-pipeline`.

**3. For ENCODE compliance, set `atac.pipeline_type` to `tf` (transcription factor) or `histone` based on the experiment type.**

## Quick Start: PEPATAC

Use PEPATAC for rapid prototyping or when a lightweight, Python-based workflow is preferred.

**1. Define the PEP project configuration (PEP format):**

```yaml
pep_version: "2.0.0"
sample_table:
  path: samples.csv
```

**2. Create the samples CSV:**

```csv
sample_name,fastq1,fastq2,genome
S1,/data/S1_R1.fastq.gz,/data/S1_R2.fastq.gz,hg38
```

**3. Run the pipeline:**

```bash
looper run pepatac_config.yaml
```

PEPATAC uses serial mitochondrial read alignment (maps reads to chrM first, then filters), which differs from the post-alignment filtering approach used by nf-core and ENCODE.

## QC Gate Thresholds (质控阈值)

After running any pipeline, evaluate sample quality against these thresholds. Flag or discard samples that fall below the "acceptable" column.

| Metric | Excellent | Acceptable | Below Threshold |
|--------|-----------|------------|-----------------|
| **TSS enrichment** | > 10 | > 6 | < 6 (poor signal) |
| **FRiP (Fraction of Reads in Peaks)** | > 0.30 | > 0.20 | < 0.20 |
| **NSC (Normalized Strand Cross-correlation)** | > 1.10 | > 1.05 | < 1.05 |
| **RSC (Relative Strand Cross-correlation)** | > 1.00 | > 0.80 | < 0.80 |
| **NRF (Nucleosome-free fraction of reads)** | > 0.90 | > 0.70 | < 0.70 |
| **PBC1 (PCR Bottlenecking Coefficient)** | > 0.70 | > 0.50 | < 0.50 |
| **Mitochondrial reads (%)** | < 5% | < 20% | > 20% (severe contamination) |
| **Duplicate rate (%)** | < 20% | < 50% | > 50% (low library complexity) |
| **IDR reproducible peaks** | > 100k | > 70k | < 70k (paired replicates) |
| **Total reads (post-filter)** | > 50M | > 25M | < 25M |

**Key interpretation notes:**

- TSS enrichment is the single most informative metric. Values below 6 indicate poor chromatin accessibility signal or failed transposition.
- FRiP correlates strongly with signal-to-noise ratio. Low FRiP often reflects high mitochondrial contamination or poor nuclei isolation.
- NSC/RSC are computed by cross-correlation analysis (phantompeakqualtools). RSC < 0.8 is a reliable indicator of poor library quality.
- PBC1 measures library complexity. Values below 0.5 suggest severe PCR duplication, indicating insufficient starting material.

## Output Overview (输出文件说明)

Each pipeline produces a standard set of output files. The table below maps expected outputs across all three pipelines.

| Output | Format | Description |
|--------|--------|-------------|
| Aligned reads | BAM / CRAM | Genome-aligned, duplicate-marked, filtered |
| Peaks | narrowPeak / broadPeak | MACS2-called peaks with q-values and signal scores |
| Consensus peaks | narrowPeak | IDR-merged reproducible peaks (nf-core, ENCODE) |
| Coverage tracks | bigWig | Normalized signal tracks for genome browser visualization |
| QC reports | HTML / JSON | MultiQC or ENCODE-format quality reports |
| Alignment stats | TXT / TSV | Mapping rate, duplication rate, insert size distribution |
| Peak QC | TSV / HTML | Fraction of reads in peaks, TSS enrichment, cross-correlation |

Typical output directory layout for nf-core/atacseq:

```
results_atac/
├── bam/
├── bowtie2/
├── macs2/
│   ├── narrow_peak/
│   └── broad_peak/
├── bigwig/
├── multiqc/
└── pipeline_info/
```

## Cross-reference to Downstream Analysis

After upstream processing is complete and samples pass QC gates, proceed to downstream analysis with the **atac-chromosome** skill. That skill covers:

- Genome browser track preparation (IGV sessions)
- Chromosome-level signal visualization
- Differential accessibility analysis (DESeq2, edgeR)
- Peak annotation and gene assignment
- Motif enrichment analysis

## Additional Resources (参考资料)

### Reference Files

Detailed documentation is available in the `references/` directory. Consult these for in-depth guidance:

- **`references/nfcore-atacseq.md`** - Complete nf-core/atacseq guide: installation, all parameters, workflow steps, output structure, genome support, and common run commands
- **`references/encode-atacseq.md`** - Complete ENCODE-DCC pipeline guide: Caper setup, JSON input format, genome TSV files, ENCODE QC gates, backend configuration, and IDR analysis
- **`references/troubleshooting.md`** - Common issues organized by stage (environment, input, alignment, QC, peak calling), recovery playbooks, and resource estimation for different experiment sizes

### Example Files

Working configuration examples in `examples/`:

- **`examples/samplesheet.csv`** - nf-core/atacseq samplesheet template with multi-replicate and multi-lane examples
- **`examples/params.yaml`** - Nextflow custom parameters file with commented sections for aligner choice, peak caller settings, and skip options
- **`examples/encode_input.json`** - ENCODE-DCC Caper input JSON for paired-end ATAC-seq with hg38 genome

### Utility Scripts

Helper scripts in `scripts/`:

- **`scripts/setup_upstream_env.sh`** - Recreate the `upstream` mamba environment from scratch with all required tools (supports `--dry-run` and `-y` flags)
- **`scripts/qc_gate_check.py`** - Parse ATAC-seq QC metrics and evaluate against ENCODE thresholds with PASS/WARN/FAIL output (supports `--json`, `--test`, and custom thresholds)
- **`scripts/validate_inputs.sh`** - Validate nf-core samplesheets and ENCODE JSON inputs before running pipelines (supports `--type nfcore|encode` and `--dry-run`)
