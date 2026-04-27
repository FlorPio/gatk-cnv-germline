# nf-core/germlinecnv

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

## Introduction

**nf-core/germlinecnv** is a bioinformatics pipeline for germline copy number variation (CNV) detection using GATK GermlineCNVCaller.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It follows [nf-core](https://nf-co.re) best practices and guidelines.

## Pipeline summary

![Pipeline Diagram](gCNV.drawio.png)

1. **BedToIntervalList** - Convert BED file to interval list format
2. **CollectReadCounts** - Collect read counts for each sample
3. **DetermineGermlineContigPloidy** - Determine contig ploidy for all samples
4. **AnnotateIntervals** - Annotate intervals with GC content and mappability
5. **GermlineCNVCaller (COHORT)** - Generate Panel of Normals
6. **GermlineCNVCaller (CASE)** - Call CNVs per sample using PON
7. **PostprocessGermlineCNVCalls** - Generate VCFs and copy ratio files
8. **Annotate CNVs** - Optional annotation with MANE transcripts

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=23.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility

3. Download the pipeline and test it on a minimal dataset:

   ```bash
   nextflow run main.nf -profile test,docker --outdir results
   ```

4. Start running your own analysis:

   ```bash
   nextflow run main.nf \
       --input samplesheet.csv \
       --fasta /path/to/hg38.fa \
       --bed /path/to/targets.bed \
       --outdir results \
       -profile docker
   ```

## Documentation

### Input samplesheet

The pipeline requires a CSV samplesheet with the following format:

```csv
sample,bam,bai
sample1,/path/to/sample1.bam,/path/to/sample1.bam.bai
sample2,/path/to/sample2.bam,/path/to/sample2.bam.bai
```

| Column | Description |
|--------|-------------|
| `sample` | Sample identifier (must be unique) |
| `bam` | Full path to BAM file |
| `bai` | Full path to BAM index (optional, will look for .bam.bai if not provided) |

### Pipeline modes

The pipeline supports three execution modes:

#### 1. Full mode (default)

Generates a Panel of Normals and then runs all samples in CASE mode:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --mode full \
    --fasta /path/to/hg38.fa \
    --bed /path/to/targets.bed \
    --outdir results \
    -profile docker
```

#### 2. PON mode

Only generates the Panel of Normals (useful for creating a reusable PON with normal samples):

```bash
nextflow run main.nf \
    --input normals.csv \
    --mode pon \
    --run_case_after_pon false \
    --fasta /path/to/hg38.fa \
    --bed /path/to/targets.bed \
    --outdir pon_results \
    -profile docker
```

#### 3. CASE mode

Runs samples against an existing PON (for processing new samples):

```bash
nextflow run main.nf \
    --input new_samples.csv \
    --mode case \
    --pon_model /path/to/pon_results/pon/pon-model \
    --ploidy_model /path/to/pon_results/ploidy/ploidy-model \
    --fasta /path/to/hg38.fa \
    --bed /path/to/targets.bed \
    --outdir case_results \
    -profile docker
```

If you already have count files:

```bash
nextflow run main.nf \
    --input new_samples.csv \
    --mode case \
    --pon_model /path/to/pon-model \
    --ploidy_model /path/to/ploidy-model \
    --counts_dir /path/to/counts \
    --outdir case_results \
    -profile docker
```

### Parameters

#### Required parameters

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to samplesheet CSV |
| `--fasta` | Path to reference FASTA (with .fai and .dict) |
| `--bed` | Path to BED file with target regions |
| `--outdir` | Output directory |

#### Pipeline mode parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--mode` | `full` | Execution mode: `pon`, `case`, or `full` |
| `--run_case_after_pon` | `true` | Run CASE mode after PON generation |

#### CASE mode parameters

| Parameter | Description |
|-----------|-------------|
| `--pon_model` | Path to existing PON model directory |
| `--ploidy_model` | Path to existing ploidy model directory |
| `--counts_dir` | Path to pre-computed counts (optional) |
| `--intervals` | Path to existing interval list (optional) |

#### Annotation parameters

| Parameter | Description |
|-----------|-------------|
| `--mane_file` | Path to MANE transcript file for annotation |
| `--allosomal_contigs` | Allosomal contigs (default: `chrX`) |

#### Resource parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | 10 | Maximum CPUs per process |
| `--max_memory` | 128.GB | Maximum memory per process |
| `--max_time` | 240.h | Maximum time per process |

### Output

```
results/
├── intervals/
│   ├── intervals.interval_list
│   └── annotated_intervals.annotated.tsv
├── counts/
│   └── *.counts.hdf5
├── ploidy/
│   ├── ploidy-calls/
│   └── ploidy-model/
├── pon/
│   ├── pon-model/
│   └── pon-calls/
├── case_calls/
│   └── {sample}/
│       └── {sample}-calls/
├── results/
│   └── {sample}/
│       ├── {sample}_denoised_copy_ratios.tsv
│       ├── {sample}_intervals.vcf.gz
│       └── {sample}_segments.vcf.gz
├── annotated/
│   └── {sample}/
│       └── {sample}_annotated.tsv
└── pipeline_info/
    ├── execution_report_*.html
    ├── execution_timeline_*.html
    ├── execution_trace_*.txt
    ├── pipeline_dag_*.svg
    └── software_versions.yml
```

## Pipeline structure

```
germlinecnv/
├── main.nf                          # Entry point
├── nextflow.config                  # Main configuration
├── nextflow_schema.json             # Parameter schema
├── workflows/
│   └── germlinecnv.nf              # Main workflow
├── subworkflows/
│   └── local/
│       ├── generate_pon/           # PON generation subworkflow
│       ├── call_cnv_case/          # CASE calling subworkflow
│       └── utils_nfcore_*/         # Utility subworkflows
├── modules/
│   └── local/
│       ├── gatk4/
│       │   ├── bedtointervallist/
│       │   ├── collectreadcounts/
│       │   ├── determinegermlinecontigploidy/
│       │   ├── annotateintervals/
│       │   ├── germlinecnvcaller_cohort/
│       │   ├── germlinecnvcaller_case/
│       │   └── postprocessgermlinecnvcalls/
│       ├── annotate_cnv_vcfs/
│       └── custom/
│           └── dumpsoftwareversions/
├── conf/
│   ├── base.config                 # Base resource configuration
│   ├── modules.config              # Module-specific configuration
│   ├── test.config                 # Test profile
│   └── test_full.config            # Full test profile
├── assets/
│   ├── contig_ploidy_priors.tsv   # Default ploidy priors
│   ├── schema_input.json          # Samplesheet schema
│   └── samplesheet.csv            # Example samplesheet
├── bin/
│   └── exon_annotation.R          # Annotation script
└── docs/
    └── usage.md                    # Detailed usage documentation
```

## Typical clinical workflow

### 1. Initial PON creation (once)

Create a PON using your laboratory's normal samples:

```bash
nextflow run main.nf \
    --input normals.csv \
    --mode pon \
    --run_case_after_pon false \
    --fasta /references/hg38.fa \
    --bed /panels/trusight_cancer.bed \
    --outdir /lab/pon_trusight \
    -profile singularity
```

### 2. Process patient samples

Run new patient samples against your established PON:

```bash
nextflow run main.nf \
    --input patient_batch.csv \
    --mode case \
    --pon_model /lab/pon_trusight/pon/pon-model \
    --ploidy_model /lab/pon_trusight/ploidy/ploidy-model \
    --fasta /references/hg38.fa \
    --bed /panels/trusight_cancer.bed \
    --mane_file /annotations/MANE.GRCh38.txt \
    --outdir /results/batch_2024_01 \
    -profile singularity
```

## Troubleshooting

### Common issues

1. **Missing BAI files**: Ensure index files exist as either `sample.bam.bai` or `sample.bai`

2. **Memory errors**: Increase memory limits:
   ```bash
   nextflow run main.nf --max_memory 256.GB ...
   ```

3. **Resume failed runs**: Use the `-resume` flag:
   ```bash
   nextflow run main.nf -resume ...
   ```

4. **Docker permission issues**: The docker profile includes user emulation. If issues persist:
   ```bash
   nextflow run main.nf -profile docker --docker.runOptions '-u $(id -u):$(id -g)' ...
   ```

## Credits

nf-core/germlinecnv was originally written by the bioinformatics team.

## Citations

If you use nf-core/germlinecnv for your analysis, please cite:

- **GATK**: McKenna A, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010;20(9):1297-1303.

- **Nextflow**: Di Tommaso P, et al. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017;35(4):316-319.

## License

This pipeline is released under the MIT License.
