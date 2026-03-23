# ENCODE-DCC ATAC-seq Pipeline Complete Guide

## Overview

The ENCODE-DCC ATAC-seq pipeline is the gold-standard processing workflow for ATAC-seq and DNase-seq data, developed and maintained by the ENCODE Data Coordination Center. Unlike nf-core/atacseq, this pipeline enforces strict ENCODE quality control gates and uses Irreproducible Discovery Rate (IDR) analysis to assess replicate consistency. It produces ENCODE-compliant outputs suitable for submission to the ENCODE portal.

The pipeline is written in WDL (Workflow Description Language) and executed via Cromwell with Caper as the CLI wrapper. It supports end-to-end processing from raw FASTQ through peak calling and signal track generation, or can start from intermediate formats (BAM, filtered BAM, TAG-ALIGN).

- **License**: MIT
- **Framework**: WDL with Cromwell backend
- **Citation**: [10.5281/zenodo.156534](https://doi.org/10.5281/zenodo.156534)
- **Source**: [github.com/ENCODE-DCC/atac-seq-pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline)
- **Protocol spec**: [ENCODE ATAC-seq processing](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit)

## System Requirements

| Resource      | Minimum        | Recommended    |
|---------------|----------------|----------------|
| CPU cores     | 8              | 16+            |
| RAM           | 32 GB          | 64 GB          |
| Storage       | 100 GB         | 200+ GB        |
| Java          | 8+             | OpenJDK 17     |
| Python        | 3.8+           | 3.10+          |
| Docker        | 20.10+         | Latest stable  |
| Singularity   | 3.8+           | 4.0+           |
| Caper         | 2.0+           | Latest         |
| Cromwell      | (bundled)      | (bundled)      |

Caper bundles Cromwell, so a separate Cromwell installation is not required. For HPC environments, ensure the backend (SLURM, SGE, PBS) is configured before running.

## Installation

### Install Caper

```bash
pip install caper
```

### Initialize Backend

Choose a backend matching the compute environment. Edit the generated config file afterward.

```bash
caper init [local|slurm|sge|pbs|google]
vi ~/.caper/default.conf
```

### Clone the Pipeline

```bash
git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
cd atac-seq-pipeline
```

### (Optional) Conda Environment

The Conda method installs pipeline-specific dependencies. Use Singularity instead if Conda conflicts arise.

```bash
bash scripts/install_conda_env.sh
```

Run with `--conda` instead of `--docker` or `--singularity` after installing the Conda environments.

## Input Formats

The pipeline accepts four starting points. All use the same JSON input format but specify different fields.

### Starting from FASTQ (paired-end)

```json
{
    "atac.pipeline_type": "atac",
    "atac.title": "My ATAC-seq experiment",
    "atac.description": "Human CD4+ T cells",
    "atac.paired_end": true,
    "atac.genome_tsv": "https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v4/hg38.tsv",
    "atac.fastqs_rep1_R1": ["/data/rep1_R1.fastq.gz"],
    "atac.fastqs_rep1_R2": ["/data/rep1_R2.fastq.gz"],
    "atac.fastqs_rep2_R1": ["/data/rep2_R1.fastq.gz"],
    "atac.fastqs_rep2_R2": ["/data/rep2_R2.fastq.gz"],
    "atac.auto_detect_adapter": true
}
```

Technical replicates are specified as additional entries in the array. They are merged automatically before processing.

```json
{
    "atac.fastqs_rep1_R1": ["rep1_techrep1_R1.fq.gz", "rep1_techrep2_R1.fq.gz"],
    "atac.fastqs_rep1_R2": ["rep1_techrep1_R2.fq.gz", "rep1_techrep2_R2.fq.gz"]
}
```

### Starting from BAM

```json
{
    "atac.pipeline_type": "atac",
    "atac.paired_end": true,
    "atac.genome_tsv": "https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v4/hg38.tsv",
    "atac.bams": ["/data/rep1.bam", "/data/rep2.bam"]
}
```

### Starting from TAG-ALIGN

```json
{
    "atac.pipeline_type": "atac",
    "atac.paired_end": true,
    "atac.genome_tsv": "https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v4/hg38.tsv",
    "atac.tas": ["/data/rep1.tagAlign.gz", "/data/rep2.tagAlign.gz"]
}
```

### Mixed Input Types

Different replicates can start from different formats using indexed arrays with null placeholders:

```json
{
    "atac.paired_ends": [true, false, true],
    "atac.fastqs_rep1_R1": ["rep1_R1.fastq.gz"],
    "atac.fastqs_rep1_R2": ["rep1_R2.fastq.gz"],
    "atac.bams": [null, "rep2.bam", null],
    "atac.tas": [null, null, "rep3.tagAlign.gz"]
}
```

## Key JSON Parameters

| Parameter | Required | Description |
|-----------|----------|-------------|
| `atac.pipeline_type` | Yes | `atac` or `dnase` |
| `atac.title` | Yes | Experiment title for QC HTML report |
| `atac.description` | Yes | Experiment description for QC HTML report |
| `atac.paired_end` | Yes | Boolean, true for all paired-end replicates |
| `atac.paired_ends` | No | Per-replicate array of booleans for mixed read types |
| `atac.genome_tsv` | Yes | Path or URL to genome TSV configuration file |
| `atac.fastqs_rep*_R1` | No | Array of R1 FASTQ files per replicate |
| `atac.fastqs_rep*_R2` | No | Array of R2 FASTQ files per replicate |
| `atac.bams` | No | Array of BAM files (one per replicate) |
| `atac.nodup_bams` | No | Array of filtered/deduped BAMs |
| `atac.tas` | No | Array of TAG-ALIGN files |
| `atac.auto_detect_adapter` | No | Auto-detect Illumina/Nextera/smallRNA adapters |
| `atac.adapter` | No | Single adapter string for all FASTQs |
| `atac.adapters_rep*_R*` | No | Per-file adapter arrays |
| `atac.align_only` | No | Skip peak calling, produce filtered BAMs only |
| `atac.true_rep_only` | No | Disable pseudo-replicate generation |
| `atac.multimapping` | No | Max multimapping locations (default 4) |
| `atac.subsample_reads` | No | Subsample to N reads (0 = no subsampling) |
| `atac.read_len` | No | Array of read lengths per replicate (for BAM/TAG-ALIGN input) |
| `atac.pval_thresh` | No | MACS2 p-value threshold |
| `atac.smooth_win` | No | MACS2 smoothing window size |
| `atac.allow_multi_mapping` | No | Allow multi-mapping reads in final output |

## Genome TSV Files

The pipeline uses genome TSV files to configure all genome-specific paths (aligner index, chromosome sizes, blacklist, TSS file). Use absolute paths or HTTPS URLs.

| Genome | URL |
|--------|-----|
| hg38 | `https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v4/hg38.tsv` |
| mm10 | `https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v4/mm10.tsv` |

For custom genomes, build a genome TSV file using the provided builder script. See `docs/build_genome_database.md` in the pipeline repo for instructions.

Individual genome parameters (e.g., blacklist) can be overridden in the input JSON without modifying the TSV:

```json
{
    "atac.genome_tsv": "/data/hg38.tsv",
    "atac.blacklist": "/data/custom_blacklist.bed.gz"
}
```

## Adapter Handling

The ENCODE pipeline provides three strategies for adapter trimming. The choice depends on library preparation kit and available metadata.

### Auto-Detection

Set `atac.auto_detect_adapter` to `true`. The pipeline scans FASTQ files for known adapter sequences:
- **Illumina**: `AGATCGGAAGAGC` (most common for standard TruSeq and Nextera kits)
- **Nextera**: `CTGTCTCTTATA` (Nextera Tn5-based libraries)
- **smallRNA**: `TGGAATTCTCGG` (small RNA sequencing adapters)

Auto-detection examines the first few thousand reads of each FASTQ file and matches against these known sequences. This is the recommended approach when the library preparation kit is unknown.

### Global Adapter Override

Set `atac.adapter` to a single adapter string. This adapter is applied to all FASTQ files in the run. Auto-detection is disabled when this parameter is set.

### Per-File Adapter Specification

Define `atac.adapters_rep*_R*` arrays matching the structure of the FASTQ arrays. Each adapter string corresponds positionally to its FASTQ file. Empty strings skip trimming for that file. This is the most granular approach and is useful when technical replicates were prepared with different kits.

## Workflow Steps

1. **Adapter trimming**: Cutadapt-based trimming. Adapters can be auto-detected (Illumina, Nextera, smallRNA) or specified manually per FASTQ.
2. **BWA alignment**: Map reads to the reference genome. The pipeline uses BWA-MEM with default parameters.
3. **MAPQ filtering**: Retain only reads with MAPQ >= 30. Low-quality alignments are discarded.
4. **Mitochondrial read removal**: Discard reads mapping to the mitochondrial chromosome (chrM).
5. **Duplicate marking**: Mark PCR duplicates with Picard MarkDuplicates.
6. **TSS enrichment calculation**: Measure signal enrichment at transcription start sites using the `compute_tss_enrichment.py` script.
7. **FRiP calculation**: Fraction of reads in peaks. Computed as reads overlapping called peaks divided by total mapped reads.
8. **NSC/RSC calculation**: Normalized strand cross-correlation and relative strand cross-correlation. Computed using `phantompeakqualtools`.
9. **Peak calling**: MACS2 peak calling with configurable p-value threshold and smoothing window.
10. **Pseudo-replicate generation**: Split reads into two pseudo-replicates per biological replicate for IDR analysis.
11. **IDR analysis**: Compare true replicates and pseudo-replicates using IDR to assess reproducibility. Produces an IDR score for each peak.
12. **Final peak set**: Consensus peaks passing IDR thresholds, merged across replicates.
13. **Signal track generation**: P-value and fold-change signal tracks in bigWig format.
14. **QC HTML report**: Comprehensive HTML report with all QC metrics, IDR plots, and TSS enrichment plots.

## IDR Analysis Explained

IDR (Irreproducible Discovery Rate) is the cornerstone of ENCODE's reproducibility assessment. It quantifies the consistency of peak calls across replicates.

### How IDR Works

For each pair of replicates (true biological replicates), IDR compares the ranked peak lists and identifies peaks that are reproducible. The pipeline also generates pseudo-replicates by randomly splitting each biological replicate's reads in half, then compares:
- True replicates (rep1 vs rep2)
- Self-pseudo-replicates (rep1_pseudo1 vs rep1_pseudo2)
- Pooled pseudo-replicates (pool of all reads, split in half)

### IDR Output Interpretation

- **IDR score**: Lower values indicate higher reproducibility. The ENCODE pipeline uses an IDR threshold of 0.05 for the final peak set.
- **Conservative peak set**: Peaks reproducible in true replicate comparisons. Requires IDR >= 0.80 rescue rate.
- **Optimal peak set**: Peaks reproducible in pooled pseudo-replicate comparisons. Less stringent, typically used for exploratory analysis.
- **IDR plot**: The standard ENCODE IDR plot shows the number of peaks passing IDR threshold at various cutoffs. A steep initial decline followed by a plateau indicates good reproducibility.

### When IDR Fails

IDR analysis requires at least 2 biological replicates with sufficient peak overlap. If replicates are biologically different (e.g., different cell types mislabeled as replicates), IDR will produce very few reproducible peaks. The pipeline will report a low rescue rate and the QC report will flag the failure. In such cases, verify the experimental design before re-running.

## QC Gates (ENCODE Standard)

The ENCODE pipeline enforces specific QC thresholds. These are informational thresholds; the pipeline reports whether each metric passes but does not halt on failure.

| Metric | Threshold | Description |
|--------|-----------|-------------|
| TSS enrichment | >= 8 (pass), 4-8 (acceptable), < 4 (fail) | Signal enrichment at transcription start sites |
| FRiP | >= 0.20 (pass), 0.01-0.20 (acceptable), < 0.01 (fail) | Fraction of reads in peaks |
| NRF (Non-Redundant Fraction) | >= 0.80 (pass), 0.50-0.80 (acceptable), < 0.50 (fail) | Fraction of non-redundant mapped reads |
| PBC1 (PCR Bottlenecking 1) | >= 0.80 (pass), 0.50-0.80 (acceptable), < 0.50 (fail) | Ratio of loci with exactly one read to loci with >=1 read |
| PBC2 (PCR Bottlenecking 2) | >= 1.0 (pass) | Ratio of loci with exactly two reads to loci with exactly one read |
| NSC (Normalized Strand Cross-correlation) | >= 1.05 | Cross-correlation normalized by auto-correlation |
| RSC (Relative Strand Cross-correlation) | >= 0.80 | Cross-correlation normalized by fragment-length cross-correlation |
| IDR reproducibility | >= 0.80 for optimal, >= 2 for conservative | Peak reproducibility between replicates |
| Peak count (IDR) | >= 70,000 (conservative), >= 100,000 (optimal) | Number of peaks passing IDR threshold |
| Mitochondrial reads | < 40% | Fraction of reads mapping to mitochondria |

A sample passes ENCODE standards when all metrics meet the "pass" threshold. Samples meeting "acceptable" thresholds may proceed with caution. Samples failing thresholds require resequencing or library preparation.

## Output Files

After running the pipeline, use [Croo](https://github.com/ENCODE-DCC/croo) to organize outputs from the Cromwell working directory into a clean structure:

```bash
pip install croo
croo /path/to/cromwell/output/metadata.json
```

The organized output includes:

| File | Description |
|------|-------------|
| `filtered.bam` | Alignment BAM with duplicates removed, low MAPQ filtered, mitochondria removed |
| `filtered.tagAlign.gz` | TAG-ALIGN format (BED) of filtered reads |
| `peaks.narrowPeak` | Narrow peaks from MACS2 |
| `peaks.broadPeak` | Broad peaks from MACS2 |
| `overlap.peaks.narrowPeak` | Overlap peak set from replicates |
| `idr.peaks.narrowPeak` | IDR-filtered peak set (reproducible peaks only) |
| `signal.pval.bigWig` | P-value signal track |
| `signal.fc.bigWig` | Fold-change signal track |
| `qc.html` | Comprehensive QC report with all metrics and plots |
| `qc.json` | Machine-readable QC metrics |
| `idr_plot.png` | IDR plot showing reproducibility |

Use `qc2tsv` to compile QC metrics from multiple experiments into a spreadsheet:

```bash
pip install qc2tsv
qc2tsv /sample1/qc.json /sample2/qc.json > qc_spreadsheet.tsv
```

## Common Run Commands

### Standard Run (Docker)

```bash
caper run atac.wdl \
  -i input.json \
  --docker \
  --max-concurrent-tasks 1
```

### Standard Run (Singularity)

```bash
caper run atac.wdl \
  -i input.json \
  --singularity \
  --max-concurrent-tasks 1
```

### HPC Run (SLURM)

Submit as a leader job. Caper manages child job submission through SLURM.

```bash
caper hpc submit atac.wdl \
  -i input.json \
  --singularity \
  --leader-job-name atac_experiment
```

Check status:

```bash
caper hpc list
```

Abort all child jobs:

```bash
caper hpc abort [JOB_ID]
```

### Alignment Only (Skip Peak Calling)

Add `atac.align_only` to the input JSON:

```json
{
    "atac.pipeline_type": "atac",
    "atac.align_only": true,
    ...
}
```

### QC Only with Custom Genome

```json
{
    "atac.pipeline_type": "atac",
    "atac.title": "Custom genome ATAC-seq",
    "atac.genome_tsv": "/data/genomes/custom/v4/custom.tsv",
    "atac.align_only": true,
    ...
}
```

## Caper Backend Configuration

### Local Backend

```bash
caper init local
vi ~/.caper/default.conf
```

Set `max_concurrent_tasks` based on available CPU cores. Each task uses 4-16 cores depending on the step. For a workstation with 16 cores, set `max_concurrent_tasks` to 2-4 to avoid resource contention.

### SLURM Backend

```bash
caper init slurm
vi ~/.caper/default.conf
```

Key configuration fields in `~/.caper/default.conf`:

- `slurm_partition`: Default partition for child jobs
- `slurm_account`: Billing account (if required)
- `slurm_time_limit`: Default wall time for child jobs
- `slurm_extra_args`: Additional sbatch arguments (e.g., `--constraint=highmem`)

The alignment step is the bottleneck and consumes the most resources. On SLURM, Caper submits child jobs for each pipeline task. The leader job manages the workflow and has minimal resource requirements. Monitor child jobs with `caper hpc list` and cancel the leader job (not individual children) to stop the entire workflow.

### Google Cloud Backend

```bash
caper init google
vi ~/.caper/default.conf
```

Requires a Google Cloud project, service account key, and Google Storage bucket for input/output. Input FASTQs and genome TSV must be accessible via `gs://` URLs. The `disk_factor` parameter controls the disk size multiplier for cloud VMs. Increase this value if large genomes or deep sequencing cause disk space errors.

### Backend-Specific Resource Overrides

Resources in the ENCODE pipeline are defined per biological replicate. The key tunable parameters are:
- `atac.align_cpu`: CPU cores for the alignment step (default: 4)
- `atac.align_mem_gb`: RAM in GB for alignment (default: 16)
- `atac.align_disk_factor`: Disk multiplier relative to input size (default: 1.5)

Override these in the input JSON only when the defaults cause resource errors. For large genomes or deep sequencing, increase `atac.align_mem_gb` to 32 or `atac.align_disk_factor` to 3.0.

## SRA Data Compatibility

Read names in paired-end FASTQ files downloaded from SRA must be consistent between R1 and R2. When using `fastq-dump` or `fasterq-dump`, do not use the `--readids` flag, which appends `.1` and `.2` suffixes to read names. Inconsistent read names cause the BWA alignment step to produce empty BAM files in the filter step.

Correct SRA download:
```bash
fasterq-dump --split-files SRR1234567
# Produces SRR1234567_1.fastq and SRR1234567_2.fastq with matching read names
```

Incorrect (produces mismatched names):
```bash
fastq-dump --readids --split-files SRR1234567
# Produces SRR1234567_1.fastq and SRR1234567_2.fastq with .1/.2 appended to read names
```
