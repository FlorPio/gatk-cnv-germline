#!/usr/bin/env Rscript
#===============================================================================
# GERMLINE CNV EXON ANNOTATION SCRIPT (v2: BED-aware)
#===============================================================================
# Annotates germline CNV VCF files with MANE Select exon information,
# restricted to exons actually present in the capture BED, and optionally 
# filters by a list of genes of interest.
#
# Input:  VCF file(s) from GATK PostprocessGermlineCNVCalls
# Output: TSV file with VCF info + gene/exon annotations
#===============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(GenomicRanges)
})

#===============================================================================
# COMMAND LINE ARGUMENTS
#===============================================================================
option_list <- list(
  make_option("--vcf",
              help = "VCF file, directory containing VCFs, or text file with list of VCF paths",
              dest = 'vcf'),
  make_option("--MANE",
              help = 'MANE Select annotation file (TSV: chr, start, end, gene_symbol, exon, transcript_id)',
              dest = 'MANE'),
  make_option("--bed",
              help = 'Capture panel BED file. MANE exons NOT overlapping this BED will be dropped from annotation.',
              dest = 'bed'),
  make_option("--min-overlap",
              help = 'Minimum bp of a MANE exon that must be covered by the BED to keep it (default: 1)',
              dest = 'min_overlap', default = 1, type = "integer"),
  make_option("--genes",
              help = 'Optional: file with list of genes to filter (one per line). If missing, NO filtering is applied.',
              dest = 'genes', default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))

vcf_file   <- opt$vcf
MANE_file  <- opt$MANE
bed_file   <- opt$bed
genes_file <- opt$genes
min_ov     <- opt$min_overlap

#===============================================================================
# VALIDATE INPUTS
#===============================================================================
if (is.null(vcf_file)  || !file.exists(vcf_file))  stop("VCF file does not exist: ",  vcf_file)
if (is.null(MANE_file) || !file.exists(MANE_file)) stop("MANE file does not exist: ", MANE_file)

use_bed <- !is.null(bed_file) && nzchar(bed_file) && file.exists(bed_file)
if (!is.null(bed_file) && !use_bed) {
  warning("BED file argument provided but file not found: ", bed_file,
          " -- continuing WITHOUT BED restriction.")
}

#===============================================================================
# DETERMINE INPUT TYPE AND GET FILE LIST
#===============================================================================
if (dir.exists(vcf_file)) {
  filenames <- list.files(vcf_file, pattern = "\\.vcf(\\.gz)?$", full.names = TRUE)
  message("Processing directory with ", length(filenames), " VCF files")
} else if (grepl("\\.vcf(\\.gz)?$", vcf_file, ignore.case = TRUE)) {
  filenames <- vcf_file
  message("Processing single VCF file")
} else {
  filenames <- trimws(readLines(vcf_file))
  filenames <- filenames[filenames != ""]
  message("Processing list of ", length(filenames), " VCF files")
}
names <- basename(filenames)

#===============================================================================
# LOAD MANE AND BED, THEN RESTRICT MANE TO EXONS COVERED BY THE BED
#===============================================================================
message("Loading MANE: ", MANE_file)
MANE_data <- read.delim(MANE_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
MANE_G_full <- makeGRangesFromDataFrame(MANE_data,
                                        seqnames.field = "chr",
                                        start.field = "start",
                                        end.field = "end",
                                        keep.extra.columns = TRUE)

if (use_bed) {
  message("Loading capture BED: ", bed_file)
  bed_df <- read.delim(bed_file, sep = "\t", header = FALSE,
                       comment.char = "#", stringsAsFactors = FALSE)
  bed_df <- bed_df[, 1:3]; colnames(bed_df) <- c("chr", "start", "end")
  # BED is 0-based half-open: convert to 1-based inclusive for GRanges
  bed_df$start <- bed_df$start + 1
  bed_G <- makeGRangesFromDataFrame(bed_df)

  # Calculate overlapping bp of each MANE exon with the BED
  ov <- findOverlaps(MANE_G_full, bed_G)
  if (length(ov) == 0) stop("No MANE exon overlaps with the BED. Check chromosome nomenclature (chr1 vs 1).")

  inter  <- pintersect(MANE_G_full[queryHits(ov)], bed_G[subjectHits(ov)])
  bp_by_exon <- tapply(width(inter), queryHits(ov), sum)

  keep_idx <- as.integer(names(bp_by_exon)[bp_by_exon >= min_ov])
  MANE_G   <- MANE_G_full[keep_idx]
  MANE_kept <- MANE_data[keep_idx, ]

  n_total   <- length(MANE_G_full)
  n_kept    <- length(MANE_G)
  n_dropped <- n_total - n_kept
  message(sprintf("MANE exons: %d total, %d kept (in BED >=%d bp), %d dropped",
                  n_total, n_kept, min_ov, n_dropped))

  dropped_df <- MANE_data[-keep_idx, , drop = FALSE]
} else {
  message("No BED restriction applied: using all MANE exons.")
  MANE_G    <- MANE_G_full
  MANE_kept <- MANE_data
  dropped_df <- MANE_data[0, , drop = FALSE]
}

#===============================================================================
# LOAD GENE LIST
#===============================================================================
genes_list <- NULL

if (!is.null(genes_file) && file.exists(genes_file)) {
  genes_list <- trimws(readLines(genes_file))
  genes_list <- genes_list[genes_list != "" & !grepl("^#", genes_list)]
  message("Loaded ", length(genes_list), " genes for filtering.")
} else {
  message("No gene list provided. No filtering will be applied.")
}

# Report which exons of genes of interest fell outside the BED (useful QC)
if (!is.null(genes_list)) {
  dropped_in_panel <- dropped_df[dropped_df$gene_symbol %in% genes_list, ]
  if (nrow(dropped_in_panel) > 0) {
    message("Exons of filtered genes outside the BED (will not be annotated):")
    tmp <- dropped_in_panel %>% arrange(gene_symbol, as.numeric(exon))
    for (k in seq_len(nrow(tmp))) {
      message(sprintf("  %s exon %s  %s:%d-%d",
                      tmp$gene_symbol[k], tmp$exon[k],
                      tmp$chr[k], tmp$start[k], tmp$end[k]))
    }
  }
}

#===============================================================================
# PROCESS EACH VCF FILE
#===============================================================================
for (i in seq_along(filenames)) {
  message("\n--- Processing: ", filenames[i], " ---")

  vcf <- read.delim(filenames[i], sep = "\t", header = FALSE,
                    comment.char = "#", stringsAsFactors = FALSE)
  colnames(vcf) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
  vcf$END <- as.integer(sub(".*END=([0-9]+).*", "\\1", vcf$INFO))
  vcf$row_id <- seq_len(nrow(vcf))

  vcf <- vcf %>%
    separate(SAMPLE,
             into = c("GT","CN","NP","QA","QS","QSE","QSS"),
             sep = ":", fill = "right")

  vcf_G <- makeGRangesFromDataFrame(vcf,
                                    seqnames.field = "CHROM",
                                    start.field = "POS", end.field = "END",
                                    keep.extra.columns = TRUE)

  # Overlap against MANE already filtered by BED
  res <- findOverlaps(vcf_G, MANE_G)
  if (length(res) == 0) {
    message("  -> No overlaps with MANE (in BED). Proceeding without annotations.")
    ann <- data.frame(row_id = integer(), gen = character(), exon = character(), stringsAsFactors = FALSE)
  } else {
    ann <- data.frame(
      row_id = vcf$row_id[queryHits(res)],
      gen    = MANE_kept$gene_symbol[subjectHits(res)],
      exon   = MANE_kept$exon[subjectHits(res)],
      stringsAsFactors = FALSE
    ) %>%
      # If a segment touches multiple exons of the same gene, collapse them into one row
      group_by(row_id, gen) %>%
      summarise(exon = paste(sort(unique(as.character(exon))), collapse = ","),
                .groups = "drop")
  }

  final <- vcf %>%
    left_join(ann, by = "row_id") %>%
    dplyr::select(CHROM, POS, END, ID, REF, ALT, QUAL, FILTER,
                  GT, CN, NP, QA, QS, QSE, QSS, gen, exon)
                  
  # Apply filter ONLY if a gene list was loaded
  if (!is.null(genes_list)) {
    final <- final %>% dplyr::filter(gen %in% genes_list)
  }

  output_file <- sub("\\.vcf.*", "_annotated.txt", filenames[i])
  write.table(final, file = output_file, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = TRUE)
  message("  -> Saved: ", output_file, "  (", nrow(final), " annotated CNVs)")
}

message("\n=== ANNOTATION COMPLETED ===")
