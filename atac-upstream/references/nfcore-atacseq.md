# nf-core/atacseq Complete Guide

## Overview

nf-core/atacseq is a community-maintained bioinformatics pipeline for ATAC-seq peak calling, annotation, and differential accessibility analysis. Built on Nextflow DSL2, it provides containerized, reproducible analysis across local workstations, HPC clusters, and cloud platforms. The pipeline covers the full upstream workflow from raw FASTQ files through quality control, alignment, post-alignment filtering, peak calling, consensus peak generation, annotation, and differential binding analysis with DESeq2.

- **License**: MIT
- **Framework**: Nextflow DSL2 (requires >=23.04.0)
- **Citation**: [10.5281/zenodo.2634132](https://doi.org/10.5281/zenodo.2634132)
- **Source**: [github.com/nf-core/atacseq](https://github.com/nf-core/atacseq)

## System Requirements

| Resource      | Minimum        | Recommended    |
|---------------|----------------|----------------|
| CPU cores     | 8              | 16             |
| RAM           | 32 GB          | 64 GB          |
| Storage       | 100 GB         | 200+ GB        |
| Docker        | 20.10+         | Latest stable  |
| Singularity   | 3.8+           | 4.0+           |
| Conda         | 4.10+          | mamba 1.0+     |
| Nextflow      | 23.04.0+       | Latest stable  |

Disk requirements scale with read count. Budget roughly 2-3x the combined FASTQ size for intermediate files and outputs. Peak calling and bigWig generation are the most memory-intensive steps; each process typically needs 8-16 GB RAM depending on genome size.

## Installation

Pull the pipeline:

```bash
nextflow pull nf-core/atacseq
```

### Docker Profile

Docker is the simplest container option for single-machine runs. Ensure the Docker daemon is running, then add `-profile docker` to any run command.

```bash
nextflow run nf-core/atacseq \
  -profile docker \
  --input samplesheet.csv \
  --outdir ./results \
  --genome GRCh38 \
  --read_length 150
```

### Singularity Profile

Singularity is preferred on HPC systems without root access. Images are automatically pulled and cached from Docker Hub.

```bash
nextflow run nf-core/atacseq \
  -profile singularity \
  --input samplesheet.csv \
  --outdir ./results \
  --genome GRCh38 \
  --read_length 150
```

To use a local Singularity image cache, set `NXF_SINGULARITY_CACHEDIR`:

```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/singularity_cache
```

### Conda Profile

Conda environments are created per process. This is slower on first run due to environment solving but avoids container overhead.

```bash
nextflow run nf-core/atacseq \
  -profile conda \
  --input samplesheet.csv \
  --outdir ./results \
  --genome GRCh38 \
  --read_length 150
```

## Input Preparation

The pipeline requires a comma-separated samplesheet (CSV) with a header row. The basic format uses four columns:

```csv
sample,fastq_1,fastq_2,replicate
CONTROL,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,1
CONTROL,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,2
TREATMENT,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,1
```

- **sample**: Condition name (used for grouping in differential analysis)
- **fastq_1**: Path to R1 FASTQ (gzip-compressed)
- **fastq_2**: Path to R2 FASTQ (leave empty for single-end)
- **replicate**: Biological replicate number within the condition

Optional columns include `adapter_type` (auto/illumina/nextera/truseq/smallrna), `trim_front1`, and `trim_front2` for specifying the number of bases to trim from read starts.

When using controls (e.g., naked DNA or input controls), add `--with_control` to the run command and list controls in the samplesheet with a distinct sample group name.

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | (required) | Path to samplesheet CSV |
| `--outdir` | (required) | Output directory |
| `--genome` | (none) | iGenomes reference name (e.g., GRCh38, GRCm38) |
| `--fasta` | (none) | Custom FASTA genome file |
| `--gtf` | (none) | Custom GTF annotation file |
| `--aligner` | `bwa` | Aligner: bwa, bowtie2, chromap, star |
| `--read_length` | (none) | Read length (50, 75, 100, 150, 200) |
| `--macs_gsize` | (auto) | Effective genome size for MACS2 |
| `--narrow_peak` | false | Call narrow peaks instead of broad |
| `--broad_cutoff` | 0.1 | Broad peak FDR cutoff |
| `--macs_fdr` | 0.05 | FDR cutoff for peak detection |
| `--macs_pvalue` | (none) | p-value cutoff (mutually exclusive with macs_fdr) |
| `--min_reps_consensus` | 1 | Min replicates for consensus peak inclusion |
| `--keep_dups` | false | Retain duplicate reads in filtered BAM |
| `--keep_multi_map` | false | Retain multi-mapping reads |
| `--keep_mito` | false | Retain mitochondrial reads |
| `--mito_name` | (auto) | Mitochondrial chromosome identifier (e.g., chrM) |
| `--blacklist` | (auto) | Path to blacklist BED file |
| `--with_control` | false | Indicates controls are in the samplesheet |
| `--skip_peak_annotation` | false | Skip HOMER peak annotation |
| `--skip_consensus_peaks` | false | Skip consensus peak generation |
| `--skip_deseq2_qc` | false | Skip DESeq2 PCA and heatmap |
| `--save_reference` | false | Save aligner index to output directory |
| `--save_trimmed` | false | Save trimmed FASTQ files |
| `--save_align_intermeds` | false | Save intermediate alignment BAMs |

## Aligner Options

The pipeline supports four aligners, selectable with `--aligner`:

### BWA (default)

BWA-MEM is the default aligner and works well for both single-end and paired-end ATAC-seq data. It provides fast alignment with good sensitivity. Adjust the minimum alignment score with `--bwa_min_score` to filter low-quality mappings (default: unset, uses BWA's internal heuristic).

### Bowtie2

Bowtie2 offers configurable sensitivity/speed tradeoffs. It supports local alignment mode which can improve mapping rates for reads with adapter contamination at read ends. Suitable for large genomes where BWA-MEM may be slow.

### Chromap

Chromap is a fast aligner optimized for chromatin accessibility data. It can process data 10-20x faster than BWA or Bowtie2. However, Chromap currently only supports paired-end reads up to the mapping step in this pipeline. Avoid using Chromap for single-end data or when downstream steps beyond alignment are needed.

### STAR

STAR is a splice-aware aligner typically used for RNA-seq but also supports ATAC-seq. Its primary advantage is handling reads spanning splice junctions, which is generally not relevant for ATAC-seq. Use STAR only if it matches an established institutional workflow or if STAR indices are already available.

## Workflow Steps

The pipeline executes the following stages in order:

1. **FastQC**: Raw read quality assessment on input FASTQ files
2. **Trim Galore**: Adapter trimming and quality filtering (wraps Cutadapt and FastQC)
3. **Alignment**: Read mapping with the selected aligner (BWA, Bowtie2, Chromap, or STAR)
4. **Picard MarkDuplicates**: Identification and marking of PCR duplicate reads
5. **Merge alignments**: Combine multiple library lanes for the same sample (Picard MergeSamFiles)
6. **Post-alignment filtering**: Remove mitochondrial reads, blacklisted regions, duplicates, secondary alignments, unmapped reads, multi-mapped reads, soft-clipped reads, reads with >4 mismatches, and (for paired-end) discordant pairs and insert sizes >2kb
7. **Picard CollectMultipleMetrics**: Alignment statistics and insert size metrics
8. **Preseq**: Library complexity estimation
9. **bigWig generation**: Normalized signal tracks scaled to 1 million mapped reads
10. **deepTools plotProfile**: Gene-body meta-profile visualization
11. **deepTools plotFingerprint**: Genome-wide enrichment estimation
12. **MACS2 peak calling**: Broad or narrow peak detection per sample and per merged replicate
13. **Consensus peaks**: BEDTools merge across all samples, with configurable reproducibility threshold
14. **featureCounts**: Read counting over consensus peaks
15. **DESeq2**: Differential accessibility analysis, PCA, and clustering
16. **HOMER**: Peak annotation relative to gene features (promoter, intron, intergenic, etc.)
17. **ataqv**: ATAC-seq specific QC report (TSS enrichment, fragment size distribution)
18. **MultiQC**: Aggregated QC report across all samples and tools
19. **IGV session**: Session file for manual visualization of bigWig tracks and peaks

## Post-Alignment Filtering Details

The filtering step is critical for ATAC-seq data quality. The pipeline applies multiple filters sequentially using SAMtools and BAMTools:

- **Mitochondrial reads**: Removed by default (controlled by `--keep_mito`). The mitochondrial chromosome name is auto-detected from iGenomes config or set manually with `--mito_name`.
- **Blacklist regions**: Reads overlapping ENCODE blacklist regions are removed. These are genomic regions with anomalous signal in sequencing assays.
- **Duplicate reads**: Removed by default (controlled by `--keep_dups`). Duplicate marking happens before filtering; removal happens during the filter step.
- **Multi-mapping reads**: Reads with MAPQ below the aligner threshold are removed (controlled by `--keep_multi_map`).
- **Secondary/supplementary alignments**: Only primary alignments are retained.
- **Unmapped reads**: Reads without a mapped position are discarded.
- **Soft-clipped reads**: Reads with soft-clipping are filtered to remove poorly aligned fragments.
- **Mismatch filter**: Reads with more than 4 mismatches are removed.
- **Insert size filter** (paired-end only): Read pairs with insert size > 2000bp are removed.
- **Chimeric filter** (paired-end only): Read pairs mapping to different chromosomes are removed.
- **Orientation filter** (paired-end only): Only FR (forward-reverse) oriented pairs are retained.
- **Concordant pair filter** (paired-end only): If one read of a pair fails any filter, the entire pair is removed.

These filters can be customized by editing the BAMTools JSON config files (`--bamtools_filter_pe_config` and `--bamtools_filter_se_config`) or by toggling individual keep flags.

## Differential Accessibility Analysis

The differential analysis module uses DESeq2 to identify regions of differential chromatin accessibility between sample conditions. This analysis requires at least two conditions with biological replicates.

Key considerations:
- The consensus peak set serves as the feature matrix. Each row is a peak, each column is a sample.
- Read counts are generated by featureCounts using the filtered BAM files.
- The DESeq2 workflow performs variance-stabilizing transformation (VST by default, or rlog if `--deseq2_vst false`), PCA, hierarchical clustering, and pairwise differential testing.
- Results include log2 fold changes, adjusted p-values (Benjamini-Hochberg), and normalized counts.
- Output plots include PCA scatter plots, sample distance heatmaps, and volcano plots for each pairwise comparison.

The analysis expects the samplesheet `sample` column to encode condition names. Replicates within the same condition are compared against replicates from other conditions. Control samples (when `--with_control` is set) are used as background for MACS2 peak calling but do not participate in the differential analysis as a separate condition.

## Output Directory Structure

```
results/
├── fastqc/
│   ├── sample1_fastqc.html
│   └── sample1_fastqc.zip
├── trim/
│   ├── sample1_val_1.fq.gz
│   └── sample1_val_2.fq.gz
├── trim_fastqc/
├── alignment/
│   ├── sample1.mapped.bam
│   └── sample1.mapped.bam.bai
├── picard_metrics/
│   └── sample1_picard_collectMultipleMetrics.txt
├── preseq/
│   └── sample1_preseq_c_curve.txt
├── bamtools_filter_stats/
├── filtered_bam/
│   ├── sample1.filtered.bam
│   └── sample1.filtered.bam.bai
├── bigwig/
│   ├── sample1.bw
│   └── sample1.filtered.bw
├── plot_profile/
├── plot_fingerprint/
├── peak_calling/
│   ├── macs2/
│   │   ├── sample1_peaks.narrowPeak
│   │   └── sample1_peaks.broadPeak
│   └── consensus_peak_set/
│       ├── consensus_peaks.bed
│       └── counts.tsv
├── annotation/
│   └── homer/
│       ├── sample1_annotated_peaks.txt
│       └── consensus_annotated_peaks.txt
├── differential_accessibility/
│   ├── deseq2/
│   │   ├── PCA.png
│   │   ├── heatmap.png
│   │   └── results.tsv
│   └── igv/
│       └── session.xml
├── ataqv/
│   └── sample1.html
└── multiqc/
    ├── multiqc_report.html
    └── multiqc_data/
```

## QC Metrics Interpretation

### MultiQC Report

The MultiQC report aggregates metrics from FastQC, Trim Galore, Picard, SAMtools, and MACS2. Key sections to inspect:

- **FastQC**: Per-base quality scores should be >30 across most of the read length. Adapter content should drop to near zero after trimming.
- **Trim Galore**: Check the percentage of reads trimmed and the post-trim length distribution.
- **Picard CollectMultipleMetrics**: Insert size distribution should show a nucleosomal ladder pattern (~200bp, ~400bp, ~600bp) characteristic of ATAC-seq.
- **Alignment rate**: Expect >70% mapping rate for good ATAC-seq libraries. Below 50% suggests contamination, wrong genome, or poor library quality.

### ataqv Report

ataqv generates ATAC-seq specific QC including:

- **TSS enrichment**: Signal enrichment at transcription start sites. Good ATAC-seq libraries show TSS enrichment >8. Below 4 indicates poor chromatin accessibility capture.
- **Fragment size distribution**: Should show clear nucleosomal periodicity (mononucleosome peak around 200bp, dinucleosome around 400bp). A dominant sub-nucleosomal peak (<100bp) is expected and indicates open chromatin regions.
- **Duplicate rate**: High duplicate rates (>50%) suggest PCR overamplification or low library complexity.
- **Mitochondrial read fraction**: ATAC-seq naturally has some mitochondrial reads due to the Tn5 transposase. Values >30% indicate over-sequencing or lysis issues.

### Preseq Complexity Curve

The Preseq library complexity curve estimates the expected yield of unique molecules. A curve that plateaus early suggests the library has been exhausted; a curve still climbing indicates additional sequencing would yield new fragments.

## Common Run Commands

### Standard Run

```bash
nextflow run nf-core/atacseq \
  -profile docker \
  --input samplesheet.csv \
  --outdir ./results \
  --genome GRCh38 \
  --read_length 150
```

### Skip Peak Calling (Alignment and QC Only)

```bash
nextflow run nf-core/atacseq \
  -profile docker \
  --input samplesheet.csv \
  --outdir ./results \
  --genome GRCh38 \
  --read_length 150 \
  --skip_consensus_peaks \
  --skip_peak_annotation \
  --skip_ataqv
```

### Custom Genome

```bash
nextflow run nf-core/atacseq \
  -profile docker \
  --input samplesheet.csv \
  --outdir ./results \
  --fasta /path/to/genome.fa \
  --gtf /path/to/annotations.gtf \
  --blacklist /path/to/blacklist.bed \
  --mito_name chrM \
  --macs_gsize 2913022398 \
  --read_length 150 \
  --save_reference
```

### With Controls and Differential Binding

```bash
nextflow run nf-core/atacseq \
  -profile docker \
  --input samplesheet_with_controls.csv \
  --outdir ./results \
  --genome GRCh38 \
  --read_length 150 \
  --with_control \
  --macs_fdr 0.05
```

### Resume a Failed Run

```bash
nextflow run nf-core/atacseq \
  -profile docker \
  --input samplesheet.csv \
  --outdir ./results \
  --genome GRCh38 \
  --read_length 150 \
  -resume
```

Add `-resume` with the same run command. Nextflow caches completed steps and only re-executes failed or modified processes.

## Genome Support

### Built-in iGenomes References

The pipeline bundles configuration for the following genomes via the iGenomes system:

| Genome   | Species           | MACS2 gsize |
|----------|-------------------|-------------|
| GRCh37   | Human (hg19)      | Pre-computed |
| GRCh38   | Human (hg38)      | Pre-computed |
| GRCm38   | Mouse (mm10)      | Pre-computed |
| GRCm39   | Mouse (mm39)      | Pre-computed |
| WBcel235 | C. elegans        | Pre-computed |
| BDGP6    | D. melanogaster   | Pre-computed |
| R64-1-1  | S. cerevisiae     | Pre-computed |
| EF2      | Zebrafish         | Pre-computed |

Blacklist regions are automatically applied for GRCh37, GRCh38, GRCm38, hg19, hg38, and mm10.

### Custom Genome

To use a non-iGenomes reference, supply `--fasta`, `--gtf`, and `--macs_gsize` directly. The pipeline will build the aligner index on the first run. Add `--save_reference` to cache the index for future runs.

```bash
nextflow run nf-core/atacseq \
  -profile docker \
  --input samplesheet.csv \
  --outdir ./results \
  --fasta /data/genomes/custom.fa \
  --gtf /data/genomes/custom.gtf \
  --macs_gsize 2500000000 \
  --mito_name MT \
  --blacklist /data/genomes/custom_blacklist.bed \
  --read_length 150 \
  --save_reference
```

The mitochondrial chromosome name must match the contig names in the FASTA file (e.g., `chrM`, `MT`, or `chr mitochondrial`). Check the FASTA headers if mitochondrial filtering is unexpectedly skipped.

### iGenomes Base Override

To point to a local or custom iGenomes mirror:

```bash
nextflow run nf-core/atacseq \
  -profile docker \
  --igenomes_base /path/to/igenomes \
  --input samplesheet.csv \
  --outdir ./results \
  --genome GRCh38 \
  --read_length 150
```
