/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOWS
//
include { GENERATE_PON    } from '../subworkflows/local/generate_pon/main'
include { CALL_CNV_CASE   } from '../subworkflows/local/call_cnv_case/main'

//
// MODULES - nf-core
//
include { GATK4_BEDTOINTERVALLIST  } from '../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_COLLECTREADCOUNTS  } from '../modules/nf-core/gatk4/collectreadcounts/main'

//
// MODULES - local
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/local/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GERMLINECNV {

    take:
    ch_samplesheet // channel: [ val(meta), path(bam), path(bai) ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Prepare reference channels with meta maps (required by nf-core modules)
    // FIX: Use Channel.value() instead of Channel.of() so channels can be reused for all samples
    //
    ch_fasta = params.fasta ? Channel.value([ [ id: 'reference' ], file(params.fasta, checkIfExists: true) ]) : Channel.empty()
    ch_fai   = params.fasta ? Channel.value([ [ id: 'reference' ], file("${params.fasta}.fai", checkIfExists: true) ]) : Channel.empty()
    ch_dict  = params.fasta ? Channel.value([ [ id: 'reference' ], file(params.fasta.replaceAll(/\.fa(sta)?$/, '.dict'), checkIfExists: true) ]) : Channel.empty()

    // Also create a simple path channel for dict (used by some processes)
    ch_dict_path = ch_dict.map { meta, dict -> dict }

    //
    // Prepare BED channel
    // FIX: Use Channel.value() instead of Channel.of()
    //
    ch_bed = params.bed ? Channel.value([ [ id: 'intervals' ], file(params.bed, checkIfExists: true) ]) : Channel.empty()

    //
    // Prepare ploidy priors channel
    //
    ch_ploidy_priors = params.ploidy_priors ? 
        file(params.ploidy_priors, checkIfExists: true) : 
        file("${projectDir}/assets/contig_ploidy_priors.tsv", checkIfExists: true)

    //
    // Prepare MANE file channel for annotation (optional)
    //
    ch_mane_file = params.mane_file ? 
        Channel.fromPath(params.mane_file, checkIfExists: true).collect() : 
        Channel.empty()

    //
    // Prepare R script channel for annotation
    //
    ch_rscript = Channel.fromPath("${projectDir}/bin/exon_annotation.R", checkIfExists: true).collect()

    //
    // Prepare genes list channel for annotation filtering (optional)
    // Uses default cancer gene panel from assets if not specified
    //
    ch_genes_list = params.genes_list ? 
        Channel.fromPath(params.genes_list, checkIfExists: true).collect() : 
        Channel.fromPath("${projectDir}/assets/cancer_genes.txt", checkIfExists: true).collect()

    //
    // WORKFLOW LOGIC BY MODE
    //
    
    if (params.mode == 'pon' || params.mode == 'full') {
        //
        // SUBWORKFLOW: Generate Panel of Normals
        //
        GENERATE_PON (
            ch_samplesheet,
            ch_bed,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_ploidy_priors
        )
        ch_versions = ch_versions.mix(GENERATE_PON.out.versions)

        // Get outputs for CASE mode
        ch_model = GENERATE_PON.out.model.map { meta, model -> model }
        ch_ploidy_model = GENERATE_PON.out.ploidy_model.map { meta, model -> model }
        ch_counts = GENERATE_PON.out.counts
        ch_intervals = GENERATE_PON.out.intervals

        if (params.mode == 'full' || params.run_case_after_pon) {
            //
            // SUBWORKFLOW: Call CNVs in CASE mode
            // Uses ploidy_model to generate sample-specific ploidy calls
            //
            CALL_CNV_CASE (
                ch_counts,
                ch_intervals,
                ch_ploidy_model,
                ch_model,
                ch_dict_path,
                ch_mane_file,
                ch_rscript,
                ch_genes_list
            )
            ch_versions = ch_versions.mix(CALL_CNV_CASE.out.versions)
        }

    } else if (params.mode == 'case') {
        //
        // MODE: CASE only (using existing PON)
        //
        if (!params.pon_model || !params.ploidy_model) {
            error "ERROR: For 'case' mode, you must provide --pon_model and --ploidy_model"
        }

        ch_model = Channel.fromPath(params.pon_model, checkIfExists: true, type: 'dir').collect()
        ch_ploidy_model = Channel.fromPath(params.ploidy_model, checkIfExists: true, type: 'dir').collect()

        // Check if we need to generate counts or use existing
        if (params.counts_dir) {
            // Use existing counts - check multiple naming conventions
            ch_counts = ch_samplesheet
                .map { meta, bam, bai ->
                    // Check all possible naming patterns
                    def counts_tsv_with_counts = file("${params.counts_dir}/${meta.id}.counts.tsv", checkIfExists: false)
                    def counts_tsv_simple = file("${params.counts_dir}/${meta.id}.tsv", checkIfExists: false)
                    def counts_hdf5_with_counts = file("${params.counts_dir}/${meta.id}.counts.hdf5", checkIfExists: false)
                    def counts_hdf5_simple = file("${params.counts_dir}/${meta.id}.hdf5", checkIfExists: false)
                    
                    def counts = null
                    if (counts_tsv_with_counts.exists()) {
                        counts = counts_tsv_with_counts
                    } else if (counts_tsv_simple.exists()) {
                        counts = counts_tsv_simple
                    } else if (counts_hdf5_with_counts.exists()) {
                        counts = counts_hdf5_with_counts
                    } else if (counts_hdf5_simple.exists()) {
                        counts = counts_hdf5_simple
                    } else {
                        error "ERROR: No counts file found for sample ${meta.id} in ${params.counts_dir}"
                    }
                    return [ meta, counts ]
                }
            ch_intervals = params.intervals ? Channel.fromPath(params.intervals, checkIfExists: true).collect() : Channel.empty()

        } else {
            // Generate counts first
            if (params.intervals) {
                ch_intervals = Channel.fromPath(params.intervals, checkIfExists: true).collect()
            } else if (params.bed) {
                GATK4_BEDTOINTERVALLIST (
                    ch_bed,
                    ch_dict
                )
                ch_intervals = GATK4_BEDTOINTERVALLIST.out.interval_list
                    .map { meta, intervals -> intervals }
                    .collect()

                ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)
            } else {
                error "ERROR: For 'case' mode without counts_dir, you must provide --bed or --intervals"
            }

            // Prepare BAM input with intervals for CollectReadCounts
            ch_bam_with_intervals = ch_samplesheet
                .combine(ch_intervals)
                .map { meta, bam, bai, intervals -> [ meta, bam, bai, intervals ] }

            GATK4_COLLECTREADCOUNTS (
                ch_bam_with_intervals,
                ch_fasta,
                ch_fai,
                ch_dict
            )
            
            ch_counts = GATK4_COLLECTREADCOUNTS.out.tsv
                .mix(GATK4_COLLECTREADCOUNTS.out.hdf5)
                .filter { meta, file -> file != null }
            
            ch_versions = ch_versions.mix(GATK4_COLLECTREADCOUNTS.out.versions)
        }

        //
        // SUBWORKFLOW: Call CNVs in CASE mode
        // Now runs DetermineGermlineContigPloidy with ploidy_model for case samples
        //
        CALL_CNV_CASE (
            ch_counts,
            ch_intervals,
            ch_ploidy_model,
            ch_model,
            ch_dict_path,
            ch_mane_file,
            ch_rscript,
            ch_genes_list
        )
        ch_versions = ch_versions.mix(CALL_CNV_CASE.out.versions)

    } else {
        error "ERROR: Invalid mode '${params.mode}'. Valid modes are: 'pon', 'case', 'full'"
    }

    //
    // Collect software versions
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    emit:
    segments       = params.mode != 'pon' ? CALL_CNV_CASE.out.segments : Channel.empty()
    multiqc_report = ch_multiqc_files
    versions       = ch_versions
}
