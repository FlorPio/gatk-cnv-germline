# Exon Annotation R Script

## Description

This R script annotates germline CNV VCF files from GATK `PostprocessGermlineCNVCalls` with exon information from MANE Select transcripts and filters results by a configurable list of genes.

## Usage

```bash
Rscript exon_annotation.R \
    --vcf <vcf_input> \
    --MANE <mane_file> \
    [--genes <genes_file>]
```

## Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--vcf` | Yes | VCF file, directory containing VCF files, or text file with list of VCF paths |
| `--MANE` | Yes | MANE Select annotation file (TSV format) |
| `--genes` | No | File with list of genes to filter (one gene per line). Default: built-in cancer panel |

## Input Formats

### VCF Input (`--vcf`)

The script accepts three input modes:

1. **Single VCF file**: Path to a `.vcf` or `.vcf.gz` file
2. **Directory**: Path to a directory containing VCF files (all `*.vcf*` files will be processed)
3. **File list**: Path to a text file containing VCF file paths (one per line)

The VCF files should be the output from GATK `PostprocessGermlineCNVCalls` with the standard format:
- Must have columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE
- INFO field must contain `END=<position>`
- SAMPLE format: `GT:CN:NP:QA:QS:QSE:QSS`

### MANE File (`--MANE`)

Tab-separated file with exon coordinates. Required columns:
```
chr	start	end	gene_symbol	exon	transcript_id
chr1	69091	70008	OR4F5	1	NM_001005484
chr1	923923	924948	SAMD11	1	NM_001385641
```

### Genes File (`--genes`)

Optional plain text file with one gene symbol per line:
```
# Comments start with #
BRCA1
BRCA2
TP53
ATM
```

If not provided, uses a default panel of 90 hereditary cancer genes.

## Output Format

Tab-separated file (`*_annotated.txt`) with columns:

| Column | Description |
|--------|-------------|
| CHROM | Chromosome |
| POS | Start position |
| END | End position |
| ID | CNV identifier |
| REF | Reference allele |
| ALT | Alternate allele (DEL, DUP, or .) |
| QUAL | Quality score |
| FILTER | Filter status |
| GT | Genotype |
| CN | Copy number |
| NP | Number of points (intervals) in segment |
| QA | Quality of all points |
| QS | Quality of segment |
| QSE | Quality of segment end |
| QSS | Quality of segment start |
| gen | Gene symbol (from MANE) |
| exon | Exon number |

## Filtering Applied

1. **Gene filter**: Only CNVs overlapping genes in the gene list are kept
2. **No quality filters**: All CNVs passing the gene filter are output regardless of quality scores

## Quality Metrics Explanation

The VCF quality metrics from GATK GermlineCNVCaller:

- **QUAL**: Phred-scaled quality score for the CNV call
- **QA**: Quality of all points - confidence across all intervals in the segment
- **QS**: Quality of segment - overall segment quality
- **QSE**: Quality of segment end - confidence at segment breakpoint (end)
- **QSS**: Quality of segment start - confidence at segment breakpoint (start)

**Recommended downstream filtering** (not applied by script):
- `QS > 100` for high-confidence calls
- `CN != 2` for actual deletions/duplications
- `NP >= 3` for segments spanning multiple intervals

## Example

```bash
# With default gene panel
Rscript exon_annotation.R \
    --vcf sample_segments.vcf \
    --MANE exons_mane.txt

# With custom gene list
Rscript exon_annotation.R \
    --vcf /path/to/vcf_directory \
    --MANE exons_mane.txt \
    --genes my_genes.txt
```

## Dependencies

R packages required (included in `florpio/cnv-annotate-r:1.0` Docker image):
- optparse
- dplyr
- tidyr
- GenomicRanges (Bioconductor)

## Output Example

```
CHROM	POS	END	ID	REF	ALT	QUAL	FILTER	GT	CN	NP	QA	QS	QSE	QSS	gen	exon
chr1	6519443	12828789	CNV_chr1_6519443_12828789	N	.	3076.53	.	0/0	2	910	9	3077	7	42	KIF1B	1
chr1	6519443	12828789	CNV_chr1_6519443_12828789	N	.	3076.53	.	0/0	2	910	9	3077	7	42	KIF1B	2
chr17	43044295	43125364	CNV_chr17_43044295_43125364	N	<DEL>	521.18	.	1/1	1	15	12	521	8	15	BRCA1	1
```
