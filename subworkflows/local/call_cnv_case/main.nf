/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL CNV CASE SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calls CNVs using CASE mode with existing PON model using nf-core modules
    
    Important notes:
    - Runs DetermineGermlineContigPloidy with --model for case samples
    - Uses sample-specific ploidy calls for GermlineCNVCaller
    - GATK4_POSTPROCESSGERMLINECNVCALLS does NOT support Conda
    - Must use Docker, Singularity, or Podman
    - Optional exon annotation with MANE transcript file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GATK4_DETERMINEGERMLINECONTIGPLOIDY } from '../../../modules/nf-core/gatk4/determinegermlinecontigploidy/main'
include { GATK4_GERMLINECNVCALLER             } from '../../../modules/nf-core/gatk4/germlinecnvcaller/main'
include { GATK4_POSTPROCESSGERMLINECNVCALLS   } from '../../../modules/nf-core/gatk4/postprocessgermlinecnvcalls/main'
include { ANNOTATE_CNV_EXONS                  } from '../../../modules/local/annotate_cnv_exons/main'

workflow CALL_CNV_CASE {

    take:
    ch_counts        // channel: [ val(meta), path(counts) ] - TSV or HDF5 count files
    ch_intervals     // channel: path(intervals)
    ch_ploidy_model  // channel: path(ploidy_model) - ploidy model directory from PoN
    ch_model         // channel: path(model) - PON CNV model directory
    ch_dict          // channel: path(dict) - sequence dictionary
    ch_mane_file     // channel: path(mane_file) - optional MANE transcript file for annotation
    ch_bed_file      // channel: path(bed) - capture panel BED file (path only, no meta) for annotation
    ch_rscript       // channel: path(rscript) - optional R annotation script
    ch_genes_list    // channel: path(genes_list) - optional gene list file for filtering

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: DetermineGermlineContigPloidy in CASE mode
    // Input: tuple val(meta), path(counts), path(bed), path(exclude_beds)
    //        tuple val(meta2), path(ploidy_model)
    //        path contig_ploidy_table
    // When ploidy_model is provided, runs in CASE mode
    // NOTE: In CASE mode, intervals should NOT be provided (they're in the model)
    // Output: case-specific ploidy calls
    //
    ch_ploidy_input = ch_counts
        .map { meta, counts ->
            // Format: [ meta, counts, intervals (empty for CASE mode), exclude_beds ]
            [ meta, counts, [], [] ]
        }

    // Create ploidy model channel with meta for module compatibility.
    // .first() converts the queue channel into a value channel so the same
    // model is broadcast to every sample in ch_counts. Without this, only the
    // first sample is processed (queue is consumed once) and the pipeline
    // silently reports "completed successfully" with most samples missing.
    ch_ploidy_model_with_meta = ch_ploidy_model
        .map { model -> [ [ id: 'ploidy_model' ], model ] }
        .first()

    GATK4_DETERMINEGERMLINECONTIGPLOIDY (
        ch_ploidy_input,
        ch_ploidy_model_with_meta,
        []  // No priors table when using model
    )
    ch_versions = ch_versions.mix(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.versions.first())

    //
    // Get case-specific ploidy calls (keep meta for proper joining)
    //
    ch_case_ploidy_calls = GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls

    //
    // MODULE: GermlineCNVCaller in CASE mode
    // Input: tuple val(meta), path(tsv), path(intervals), path(ploidy), path(model)
    // When model is provided, runs in CASE mode
    // Output: casecalls directory
    //
    ch_cnv_input = ch_counts
        .join(ch_case_ploidy_calls)
        .combine(ch_intervals)
        .combine(ch_model)
        .map { meta, counts, ploidy, intervals, model ->
            // Note: annotated_intervals is empty for CASE mode - model already has this info
            [ meta, counts, intervals, ploidy, model, [] ]
        }

    GATK4_GERMLINECNVCALLER (
        ch_cnv_input
    )
    ch_versions = ch_versions.mix(GATK4_GERMLINECNVCALLER.out.versions.first())

    //
    // MODULE: Postprocess CNV calls
    // Input: tuple val(meta), path(calls), path(model), path(ploidy)
    // Output: intervals VCF, segments VCF, denoised copy ratios
    //
    // NOTE: This module does NOT support Conda!
    // NOTE: PostprocessGermlineCNVCalls uses the CASE ploidy calls, not PoN ploidy calls
    // Use join() to properly match ploidy calls by sample meta.id
    //
    ch_postprocess_input = GATK4_GERMLINECNVCALLER.out.casecalls
        .join(ch_case_ploidy_calls)
        .combine(ch_model)
        .map { meta, calls, ploidy, model ->
            [ meta, calls, model, ploidy ]
        }

    GATK4_POSTPROCESSGERMLINECNVCALLS (
        ch_postprocess_input
    )
    ch_versions = ch_versions.mix(GATK4_POSTPROCESSGERMLINECNVCALLS.out.versions.first())

    //
    // MODULE: Annotate CNV exons (optional)
    // Input: tuple val(meta), path(vcf) ; path(mane_file) ; path(rscript)
    // Output: annotated TSV with exon information
    //
    ch_annotated = Channel.empty()
    
    if (ch_mane_file) {
        ANNOTATE_CNV_EXONS (
            GATK4_POSTPROCESSGERMLINECNVCALLS.out.segments,
            ch_mane_file,
            ch_bed_file,
            ch_rscript,
            ch_genes_list
        )
        ch_annotated = ANNOTATE_CNV_EXONS.out.annotated
        ch_versions = ch_versions.mix(ANNOTATE_CNV_EXONS.out.versions.first())
    }

    emit:
    ploidy_calls  = GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls      // channel: [ val(meta), path(calls) ]
    calls         = GATK4_GERMLINECNVCALLER.out.casecalls              // channel: [ val(meta), path(calls) ]
    segments      = GATK4_POSTPROCESSGERMLINECNVCALLS.out.segments     // channel: [ val(meta), path(vcf) ]
    intervals_vcf = GATK4_POSTPROCESSGERMLINECNVCALLS.out.intervals    // channel: [ val(meta), path(vcf) ]
    denoised      = GATK4_POSTPROCESSGERMLINECNVCALLS.out.denoised     // channel: [ val(meta), path(vcf) ]
    annotated     = ch_annotated                                       // channel: [ val(meta), path(tsv) ]
    versions      = ch_versions                                        // channel: path(versions.yml)
}
