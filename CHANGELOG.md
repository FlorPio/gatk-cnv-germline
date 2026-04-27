# nf-core/germlinecnv: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - 2025-02-23

Initial release of nf-core/germlinecnv.

### Features

- Three execution modes: `pon` (Panel of Normals generation), `case` (CNV calling with existing PON), `full` (both)
- GATK4 GermlineCNVCaller pipeline with nf-core modules:
  - BedToIntervalList
  - CollectReadCounts
  - AnnotateIntervals
  - DetermineGermlineContigPloidy
  - GermlineCNVCaller (COHORT and CASE modes)
  - PostprocessGermlineCNVCalls
- Support for pre-computed read counts (`--counts_dir`)
- Optional exon-level CNV annotation with MANE transcripts
- Configurable allosomal contig handling
- Custom cancer gene panel filtering
