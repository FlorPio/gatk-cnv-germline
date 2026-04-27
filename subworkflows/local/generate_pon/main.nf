/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENERATE PON SUBWORKFLOW (con soporte para counts pre-computados)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Generates Panel of Normals for germline CNV calling using nf-core modules
    
    NUEVO: Soporta --counts_dir para reutilizar read counts existentes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GATK4_BEDTOINTERVALLIST             } from '../../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_COLLECTREADCOUNTS             } from '../../../modules/nf-core/gatk4/collectreadcounts/main'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY } from '../../../modules/nf-core/gatk4/determinegermlinecontigploidy/main'
include { GATK4_ANNOTATEINTERVALS             } from '../../../modules/nf-core/gatk4/annotateintervals/main'
include { GATK4_GERMLINECNVCALLER             } from '../../../modules/nf-core/gatk4/germlinecnvcaller/main'

workflow GENERATE_PON {

    take:
    ch_bam           // channel: [ val(meta), path(bam), path(bai) ]
    ch_bed           // channel: [ val(meta), path(bed) ]
    ch_fasta         // channel: [ val(meta), path(fasta) ]
    ch_fai           // channel: [ val(meta), path(fai) ]
    ch_dict          // channel: [ val(meta), path(dict) ]
    ch_ploidy_priors // channel: path(ploidy_priors)

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Convert BED to interval list
    //
    GATK4_BEDTOINTERVALLIST (
        ch_bed,
        ch_dict
    )
    ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions.first())

    //
    // Get intervals (remove meta for use in multiple processes)
    //
    ch_intervals = GATK4_BEDTOINTERVALLIST.out.interval_list
        .first()
        .map { meta, intervals -> intervals }

    //
    // LOGIC: Check if we should use precomputed counts or generate new ones
    //
    if (params.counts_dir) {
        //
        // USE PRECOMPUTED COUNTS
        // Load .hdf5 or .tsv files from the specified directory
        //
        log.info "=========================================="
        log.info "USANDO COUNTS PRE-COMPUTADOS"
        log.info "Directorio: ${params.counts_dir}"
        log.info "=========================================="

        ch_counts = ch_bam
            .map { meta, bam, bai ->
                // Try to find the counts file (HDF5 or TSV)
                def counts_hdf5 = file("${params.counts_dir}/${meta.id}.hdf5", checkIfExists: false)
                def counts_tsv = file("${params.counts_dir}/${meta.id}.tsv", checkIfExists: false)
                def counts_hdf5_alt = file("${params.counts_dir}/${meta.id}.counts.hdf5", checkIfExists: false)
                def counts_tsv_alt = file("${params.counts_dir}/${meta.id}.counts.tsv", checkIfExists: false)
                
                def counts = null
                if (counts_hdf5.exists()) {
                    counts = counts_hdf5
                } else if (counts_tsv.exists()) {
                    counts = counts_tsv
                } else if (counts_hdf5_alt.exists()) {
                    counts = counts_hdf5_alt
                } else if (counts_tsv_alt.exists()) {
                    counts = counts_tsv_alt
                } else {
                    error "ERROR: No se encontró archivo de counts para ${meta.id} en ${params.counts_dir}. " +
                          "Buscado: ${meta.id}.hdf5, ${meta.id}.tsv, ${meta.id}.counts.hdf5, ${meta.id}.counts.tsv"
                }
                
                log.info "  Usando counts existente: ${counts}"
                return [ meta, counts ]
            }

    } else {
        //
        // GENERATE NEW COUNTS
        // Run CollectReadCounts for each sample
        //
        log.info "Generando nuevos read counts..."

        ch_bam_with_intervals = ch_bam
            .combine(ch_intervals)
            .map { meta, bam, bai, intervals -> [ meta, bam, bai, intervals ] }

        GATK4_COLLECTREADCOUNTS (
            ch_bam_with_intervals,
            ch_fasta,
            ch_fai,
            ch_dict
        )
        ch_versions = ch_versions.mix(GATK4_COLLECTREADCOUNTS.out.versions.first())

        //
        // Collect all counts for cohort analysis
        //
        ch_counts = GATK4_COLLECTREADCOUNTS.out.tsv
            .mix(GATK4_COLLECTREADCOUNTS.out.hdf5)
            .filter { meta, file -> file != null }
    }

    //
    // Collect all counts for cohort analysis
    //
    ch_all_counts = ch_counts
        .map { meta, counts -> counts }
        .collect()
        .map { counts -> [ [ id: 'cohort' ], counts ] }

    //
    // MODULE: Determine germline contig ploidy (COHORT mode)
    //
    ch_ploidy_input = ch_all_counts
        .combine(ch_intervals)
        .map { meta, counts, intervals -> [ meta, counts, intervals, [] ] }

    ch_ploidy_model_empty = Channel.value( [ [ id: 'ploidy_model' ], [] ] )

    GATK4_DETERMINEGERMLINECONTIGPLOIDY (
        ch_ploidy_input,
        ch_ploidy_model_empty,
        ch_ploidy_priors
    )
    ch_versions = ch_versions.mix(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.versions.first())

    //
    // MODULE: Annotate intervals
    //
    ch_intervals_meta = ch_intervals.map { intervals -> [ [ id: 'intervals' ], intervals ] }
    
    ch_dummy = Channel.value( [ [ id: 'null' ], [] ] )
    
    GATK4_ANNOTATEINTERVALS (
        ch_intervals_meta,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_dummy,
        ch_dummy,
        ch_dummy,
        ch_dummy
    )

    //
    // Get ploidy calls directory
    //
    ch_ploidy_calls = GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls
        .map { meta, calls -> calls }

    ch_annotated = GATK4_ANNOTATEINTERVALS.out.annotated_intervals
        .map { meta, annotated -> annotated }

    //
    // MODULE: GermlineCNVCaller in COHORT mode
    //
    ch_cnvcaller_input = ch_all_counts
        .combine(ch_intervals)
        .combine(ch_ploidy_calls)
        .combine(ch_annotated)
        .map { meta, counts, intervals, ploidy, annotated -> [ meta, counts, intervals, ploidy, [], annotated ] }

    GATK4_GERMLINECNVCALLER (
        ch_cnvcaller_input
    )
    ch_versions = ch_versions.mix(GATK4_GERMLINECNVCALLER.out.versions.first())

    emit:
    model         = GATK4_GERMLINECNVCALLER.out.cohortmodel
    cohort_calls  = GATK4_GERMLINECNVCALLER.out.cohortcalls
    ploidy_calls  = GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls
    ploidy_model  = GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.model
    intervals     = ch_intervals
    counts        = ch_counts
    annotated     = ch_annotated
    versions      = ch_versions
}
