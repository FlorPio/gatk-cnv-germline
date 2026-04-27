# nf-core/germlinecnv: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

1. **BedToIntervalList** - Convert BED to Picard interval list
2. **CollectReadCounts** - Collect read counts per interval
3. **AnnotateIntervals** - Annotate intervals with GC content and mappability
4. **DetermineGermlineContigPloidy** - Determine contig ploidy per sample
5. **GermlineCNVCaller** - Call germline CNVs (COHORT or CASE mode)
6. **PostprocessGermlineCNVCalls** - Generate per-sample VCF files
7. **AnnotateCNVExons** - (Optional) Annotate CNVs with exon information

## Output directories

### `intervals/`

- `intervals.interval_list`: Picard-format interval list converted from BED
- `annotated_intervals.tsv`: Intervals annotated with GC content

### `counts/`

- `<sample>.counts.tsv`: Read counts per interval for each sample

### `ploidy/`

- Ploidy model and per-sample ploidy calls

### `cnv_calls/`

- CNV model (COHORT mode) or per-sample CNV calls (CASE mode)

### `results/<sample>/`

- `<sample>.segments.vcf.gz`: Segmented CNV calls
- `<sample>.intervals.vcf.gz`: Interval-level CNV calls
- `<sample>.denoised_copy_ratios.tsv`: Denoised copy ratios

### `annotated/<sample>/`

- `<sample>_annotated.txt`: CNV calls annotated with exon-level information (when `--mane_file` is provided)

### `pipeline_info/`

- `software_versions.yml`: Software versions used in the run
- `execution_report_*.html`: Nextflow execution report
- `execution_timeline_*.html`: Nextflow timeline
- `execution_trace_*.txt`: Nextflow trace file
- `pipeline_dag_*.svg`: Pipeline DAG visualization
