# nf-core/germlinecnv: Usage

## Introduction

nf-core/germlinecnv is a bioinformatics pipeline for germline copy number variant (CNV) detection using GATK4 GermlineCNVCaller.

## Samplesheet input

You will need a samplesheet with information about the samples you want to analyse. It has to be a comma-separated file with 3 columns and a header row:

```csv
sample,bam,bai
SAMPLE_1,/path/to/sample1.bam,/path/to/sample1.bam.bai
SAMPLE_2,/path/to/sample2.bam,/path/to/sample2.bam.bai
```

| Column   | Description                            |
| -------- | -------------------------------------- |
| `sample` | Custom sample name.                    |
| `bam`    | Full path to BAM file.                 |
| `bai`    | Full path to BAM index file.           |

## Pipeline modes

### Full mode (default)

Generates a Panel of Normals and then runs CNV calling on all samples:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --fasta reference.fa \
    --bed targets.bed \
    --outdir results \
    -profile docker
```

### PON mode

Generates only the Panel of Normals:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --fasta reference.fa \
    --bed targets.bed \
    --mode pon \
    --outdir results \
    -profile docker
```

### CASE mode

Runs CNV calling using an existing Panel of Normals:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --fasta reference.fa \
    --pon_model /path/to/pon_model \
    --ploidy_model /path/to/ploidy_model \
    --mode case \
    --outdir results \
    -profile docker
```

## Optional parameters

- `--counts_dir`: Path to pre-computed read counts directory (skips CollectReadCounts)
- `--mane_file`: Path to MANE transcript file for exon-level annotation
- `--genes_list`: Custom gene list for annotation filtering (default: cancer gene panel)
- `--allosomal_contigs`: Comma-separated allosomal contigs (default: `chrX`)
- `--ploidy_priors`: Custom contig ploidy priors file
