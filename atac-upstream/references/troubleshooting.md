# ATAC-seq Upstream Pipeline Troubleshooting

## Environment Issues

### Java Not Found (ENCODE pipeline)

**Symptom**: Caper or Cromwell fails to start with "java: command not found" or similar.

**Root cause**: Java 8+ is required for Cromwell. Some HPC systems do not expose Java by default.

**Resolution**:
```bash
module load java/openjdk-17
# or
export JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
export PATH=$JAVA_HOME/bin:$PATH
```

Verify with `java -version`. The ENCODE pipeline requires Java 8+ but runs with OpenJDK 17.

**Prevention**: Add the Java module load to `~/.bashrc` or the Caper config file's `java` field.

### Docker Not Running

**Symptom**: nf-core pipeline fails with "Cannot connect to the Docker daemon."

**Root cause**: The Docker service is not started or the current user lacks permissions.

**Resolution**:
```bash
sudo systemctl start docker
sudo usermod -aG docker $USER
# Log out and log back in for group change to take effect
```

**Prevention**: Enable Docker to start on boot with `sudo systemctl enable docker`.

### Conda Environment Conflicts

**Symptom**: ModuleNotFoundError, version conflicts, or missing packages when using the Conda profile.

**Root cause**: A shared Conda installation with incompatible base environment, or channel priority conflicts.

**Resolution**: Use mamba for faster solving, and always create a fresh environment:
```bash
conda create -n nf-atac -c conda-forge -c bioconda nextflow nf-core
conda activate nf-atac
```

For the ENCODE pipeline, do not use a shared Conda. Install Miniconda3 from scratch and run `bash scripts/install_conda_env.sh` on an interactive compute node (not a login node).

### Nextflow Version Mismatch

**Symptom**: Pipeline fails with "This pipeline requires Nextflow version >= 23.04.0" or DSL2 syntax errors.

**Root cause**: Installed Nextflow is older than the pipeline requires.

**Resolution**:
```bash
nextflow self-update
# or install a specific version
curl -s https://get.nextflow.io | bash
```

Check version with `nextflow -v`.

## Input Issues

### Samplesheet Format Errors (nf-core)

**Symptom**: "Missing samplesheet definition" or "Malformed samplesheet" error.

**Root cause**: Wrong delimiter (tabs instead of commas), missing header row, or incorrect column names.

**Resolution**: Ensure the file is CSV (comma-separated) with this exact header:
```
sample,fastq_1,fastq_2,replicate
```

Common mistakes: using TSV, including spaces after commas, having Windows-style line endings. Fix with:
```bash
dos2unix samplesheet.csv
file samplesheet.csv  # Should show "ASCII text"
```

### Missing FASTQ Files

**Symptom**: "No such file or directory" for FASTQ paths in the samplesheet or JSON.

**Root cause**: Relative paths instead of absolute paths, or typos in filenames.

**Resolution**: Use absolute paths in all input files. Verify existence before running:
```bash
awk -F',' 'NR>1 {print $2; print $3}' samplesheet.csv | xargs ls -la
```

### Wrong Adapter Specification (ENCODE)

**Symptom**: Reads are trimmed to zero length, or adapter content remains high after trimming.

**Root cause**: Specified adapter does not match the library preparation kit used. Nextera libraries need a different adapter than Illumina TruSeq.

**Resolution**: Enable auto-detection for unknown kits:
```json
{"atac.auto_detect_adapter": true}
```

For manual specification, the common adapters are:
- Illumina: `AGATCGGAAGAGC`
- Nextera: `CTGTCTCTTATA`
- smallRNA: `TGGAATTCTCGG`

## Alignment Issues

### Low Mapping Rate

**Symptom**: Alignment rate < 50% reported by Picard or SAMtools flagstat.

**Root cause**: Wrong reference genome, contaminated sample, or poor library quality.

**Resolution steps**:
1. Verify the genome matches the species. Check FASTQ content with `fastqc` and compare against expected k-mer profiles.
2. Check for contamination by aligning a subset of reads against a multi-species database or running `kraken2`.
3. Inspect the unmapped reads with `samtools view -f 4` to identify common contaminants (adapter dimers, rRNA, spike-ins).
4. If the rate is 50-70%, consider whether the genome build is correct (e.g., GRCh37 vs GRCh38).

### Chromosome Naming Mismatch

**Symptom**: Nearly zero alignment despite correct species, or Picard errors about sequence dictionary mismatches.

**Root cause**: The reference FASTA uses `chr1` prefixes but the aligner expects `1`, or vice versa.

**Resolution**: Check the reference FASTA header:
```bash
samtools faidx genome.fa | head -5
```

If the pipeline uses `chr`-prefixed names but the aligner index was built without prefixes (or the reverse), rebuild the index or rename chromosomes. For nf-core, this is handled automatically when using iGenomes references.

### BWA/Bowtie2 Index Errors

**Symptom**: "Can't find index" or "fread failed" during alignment.

**Root cause**: Missing or corrupted aligner index files.

**Resolution**: For nf-core, add `--save_reference` on the first run to cache the index. For ENCODE, verify the genome TSV file points to valid index paths. Rebuild if necessary:
```bash
bwa index genome.fa
# or
bowtie2-build genome.fa genome_prefix
```

## QC Failures

### Low TSS Enrichment

**Symptom**: TSS enrichment < 4 (nf-core ataqv) or < 8 (ENCODE fail threshold).

**Root cause**: Poor chromatin accessibility capture, insufficient sequencing depth, or incorrect cell type.

**Resolution**:
1. Check fragment size distribution. A proper ATAC-seq library shows nucleosomal periodicity. A flat distribution without periodicity indicates the Tn5 reaction failed.
2. Verify that the correct cell type was used and that nuclei were properly isolated.
3. Increase sequencing depth if library complexity is adequate (Preseq curve still climbing).
4. Consider re-sequencing from a new library if the nucleosomal pattern is absent.

### High Mitochondrial Reads

**Symptom**: > 30% mitochondrial reads (nf-core) or > 40% (ENCODE threshold).

**Root cause**: Over-lysis during nuclei preparation, damaged cell membranes, or excessive Tn5.

**Resolution**:
1. Reduce lysis time or detergent concentration during nuclei preparation.
2. Filter mitochondria and re-analyze the remaining data (nf-core: `--keep_mito false`, which is the default).
3. If > 70% of reads are mitochondrial, the library is unlikely to be salvageable.

### Low FRiP Score

**Symptom**: FRiP < 0.01 (ENCODE fail) or very few reads falling in called peaks.

**Root cause**: Low signal-to-noise ratio, batch effects, or incorrect peak calling parameters.

**Resolution**:
1. Check alignment quality first. Poor alignment leads to low FRiP.
2. Verify the MACS2 parameters. For ATAC-seq, the `--nomodel --shift -100 --extsize 200` settings are standard. The ENCODE pipeline applies these automatically.
3. Ensure sufficient sequencing depth. ATAC-seq typically needs 50 million paired-end reads per sample for robust peak calling.

### PCR Duplication Overload

**Symptom**: Duplicate rate > 50% and Preseq curve plateaus early.

**Root cause**: Too few unique molecules in the library, excessive PCR amplification cycles.

**Resolution**:
1. Do not increase PCR cycles to compensate for low input. Instead, start with more cells.
2. Use `--keep_dups` (nf-core) only as a diagnostic, not for final analysis.
3. Consider re-sequencing from a new library with higher input cell count.

## Peak Calling Issues

### MACS2 Memory Error

**Symptom**: "MemoryError" or killed process during MACS2 peak calling.

**Root cause**: Insufficient RAM for the genome size. Human genomes require 16-32 GB for MACS2.

**Resolution**: For nf-core, the container handles memory allocation. Increase the process memory in a custom config:
```groovy
process {
    withName: 'MACS2_*' {
        memory = '64.GB'
    }
}
```

For standalone MACS2, increase available memory or reduce the genome effective size if using a custom genome.

### No Peaks Called

**Symptom**: MACS2 produces zero or near-zero peaks.

**Root cause**: Poor alignment quality, aggressive filtering, or p-value/FDR threshold set too stringently.

**Resolution**:
1. Check that the filtered BAM contains reads: `samtools view -c filtered.bam`
2. Relax the FDR threshold: `--macs_fdr 0.1` or `--macs_pvalue 0.01`
3. Verify the effective genome size matches the actual genome build.

### Too Many Peaks

**Symptom**: Hundreds of thousands of peaks, many in repetitive or blacklisted regions.

**Root cause**: Missing or incorrect blacklist, or loose p-value threshold.

**Resolution**: Ensure a blacklist is specified. For nf-core, blacklists for common genomes are auto-applied. Verify with `--blacklist` pointing to the correct BED file. Tighten the FDR threshold to 0.01.

## nf-core Specific Issues

### Container Pull Failures

**Symptom**: "Error pulling image" or Docker Hub rate limiting.

**Root cause**: Docker Hub rate limits (100 pulls/6 hours for anonymous users) or network connectivity issues.

**Resolution**:
1. Log in to Docker Hub: `docker login`
2. Use Singularity instead, which caches images locally
3. Set a Docker mirror in daemon.json or use `containerRegistry` in Nextflow config

### Resume Failures

**Symptom**: `-resume` re-executes all processes instead of resuming from the last successful step.

**Root cause**: Changed parameters, modified input files, or cleared the Nextflow work directory.

**Resolution**: Ensure the exact same command (including all parameters) is used with `-resume`. Do not modify the samplesheet between runs. If the work directory was cleaned, resume is not possible.

### Work Directory Cleanup

**Symptom**: Disk full due to accumulated intermediate files in the Nextflow work directory.

**Resolution**: Clean up completed processes while keeping the cache for resume:
```bash
nextflow clean -f
# or remove specific run's work directory
rm -rf work/5a/1b2c3d...
```

## ENCODE Specific Issues

### Cromwell Backend Errors

**Symptom**: Caper fails to submit jobs with "Cromwell server not responding" or connection refused.

**Root cause**: Cromwell not started, wrong port, or stale lock files.

**Resolution**:
```bash
caper clean
caper run atac.wdl -i input.json --docker
```

If using a local backend, verify no stale Java/Cromwell processes are running:
```bash
ps aux | grep cromwell
kill -9 [stale_pid]
```

### IDR Failures

**Symptom**: IDR analysis produces no output or fails with "empty success array."

**Root cause**: Too few peaks from one or more replicates for meaningful IDR comparison. IDR requires at least 2 biological replicates with sufficient peak overlap.

**Resolution**:
1. Verify each replicate has at least 10,000 peaks. Check individual MACS2 outputs.
2. Relax the MACS2 p-value threshold to increase peak count.
3. If replicates are truly dissimilar (different conditions mislabeled as replicates), IDR will legitimately fail. Verify the experimental design.

### Caper Configuration Issues

**Symptom**: Jobs are submitted but fail immediately with "resource not available" or queue errors.

**Root cause**: Incorrect SLURM partition, missing account, or incompatible resource requests in `~/.caper/default.conf`.

**Resolution**: Edit `~/.caper/default.conf` to match the HPC environment. Set `slurm_partition`, `slurm_account`, and `slurm_time_limit` to valid values. Test with a simple job first:
```bash
caper run atac.wdl -i test_input.json --docker --max-concurrent-tasks 1
```

## PEPATAC Specific Issues

### Refgenie Not Initialized

**Symptom**: "Refgenie genome not found" or "genome asset missing."

**Root cause**: Refgenie is not initialized or the requested genome assets have not been pulled.

**Resolution**:
```bash
refgenie init
refgenie pull hg38/fasta
refgenie pull hg38/gene_anno
refgenie pull hg38/blacklist
```

List available assets with `refgenie list`.

### Missing Genome Assets

**Symptom**: "Asset not found for genome" error during alignment or peak annotation.

**Root cause**: The specific Refgenie asset required by PEPATAC is not downloaded.

**Resolution**: Pull all required assets for the target genome. PEPATAC typically needs: fasta, bowtie2_index, gene_anno, blacklist, and tss_annotation. Check the PEPATAC documentation for the exact asset names for each pipeline version.

### PEP Format Errors

**Symptom**: "Schema validation failed" or "Missing required columns" when running with Looper.

**Root cause**: The PEP (Portable Encodable Projects) configuration file does not match the expected schema.

**Resolution**: Validate the PEP file against the schema:
```bash
peppy validate project_config.yaml
```

Ensure required columns are present in the sample annotation CSV. PEPATAC requires at minimum: sample_name, organism, genome, and the FASTQ paths.

## Recovery Playbook

### nf-core: Resume After Failure

```bash
# Use the exact same command with -resume
nextflow run nf-core/atacseq \
  -profile docker \
  --input samplesheet.csv \
  --outdir ./results \
  --genome GRCh38 \
  --read_length 150 \
  -resume
```

If a specific process consistently fails, inspect the `.command.sh` and `.command.err` files in the work directory:
```bash
find work/ -name ".command.err" -exec grep -l "ERROR" {} \;
```

### ENCODE: Restart After Failure

```bash
# Caper automatically tracks completed tasks
caper run atac.wdl -i input.json --docker
```

Caper and Cromwell cache completed task outputs. Re-running the same command resumes from the last failure point. If Cromwell metadata is corrupted, run `caper clean` first.

### Looper (PEPATAC): Rerun Failed Samples

```bash
looper rerun project_config.yaml --failed
```

This reruns only samples that failed or were not completed. For a full clean rerun:
```bash
looper run project_config.yaml --force
```

## Resource Estimation Guide

Estimate compute resources based on total reads per sample:

| Reads per sample | CPU | RAM per process | Disk (intermediate + output) |
|------------------|-----|-----------------|------------------------------|
| 5 million | 4 | 16 GB | 20 GB |
| 50 million | 8 | 32 GB | 80 GB |
| 200 million | 16 | 64 GB | 250 GB |

MACS2 peak calling is the most memory-intensive step, scaling roughly linearly with genome size and read count. Alignment (BWA/Bowtie2) is the most CPU-intensive step and benefits from multiple threads.

For HPC job scheduling, alignment tasks typically need 4-8 cores and 32 GB RAM for human genomes at 50 million reads. Post-alignment QC and peak calling can usually run with 4 cores and 16 GB RAM. Preseq and ataqv need 2 cores and 8 GB RAM.
