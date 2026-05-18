/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENERATE PON SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    - Hybrid CollectReadCounts: reuse files from --counts_dir, generate the rest.
    - Optional FilterIntervals (params.use_filter_intervals).
    - Emits a pon_manifest.json documenting samples, BED, FASTA, intervals.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GATK4_BEDTOINTERVALLIST             } from '../../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_PREPROCESSINTERVALS           } from '../../../modules/local/gatk4/preprocessintervals/main'
include { GATK4_COLLECTREADCOUNTS             } from '../../../modules/nf-core/gatk4/collectreadcounts/main'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY } from '../../../modules/nf-core/gatk4/determinegermlinecontigploidy/main'
include { GATK4_ANNOTATEINTERVALS             } from '../../../modules/nf-core/gatk4/annotateintervals/main'
include { GATK4_FILTERINTERVALS              } from '../../../modules/nf-core/gatk4/filterintervals/main'
include { GATK4_GERMLINECNVCALLER             } from '../../../modules/nf-core/gatk4/germlinecnvcaller/main'
include { PON_MANIFEST                        } from '../../../modules/local/pon_manifest/main'

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
    // BedToIntervalList -> Picard interval_list (just a format conversion).
    // PreprocessIntervals -> applies --padding (params.bed_padding) and the
    // OVERLAPPING_ONLY merging rule. This is the canonical GATK CNV step
    // and the right place for padding (BedToIntervalList does NOT support
    // a padding argument).
    //
    GATK4_BEDTOINTERVALLIST ( ch_bed, ch_dict )
    ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions.first())

    GATK4_PREPROCESSINTERVALS (
        GATK4_BEDTOINTERVALLIST.out.interval_list,
        ch_fasta,
        ch_fai,
        ch_dict,
        Channel.value( [ [ id: 'null' ], [] ] )
    )
    ch_versions = ch_versions.mix(GATK4_PREPROCESSINTERVALS.out.versions.first())

    ch_intervals = GATK4_PREPROCESSINTERVALS.out.interval_list
        .first()
        .map { meta, intervals -> intervals }

    //
    // Helper: lookup precomputed counts file in --counts_dir
    //
    def find_existing_counts = { sample_id ->
        if (!params.counts_dir) return null
        def candidates = [
            file("${params.counts_dir}/${sample_id}.hdf5",        checkIfExists: false),
            file("${params.counts_dir}/${sample_id}.tsv",         checkIfExists: false),
            file("${params.counts_dir}/${sample_id}.counts.hdf5", checkIfExists: false),
            file("${params.counts_dir}/${sample_id}.counts.tsv",  checkIfExists: false)
        ]
        return candidates.find { it.exists() }
    }

    //
    // HYBRID counts: split samples into [existing] / [missing]
    //
    ch_bam_branched = ch_bam
        .map { meta, bam, bai -> [ meta, bam, bai, find_existing_counts(meta.id) ] }
        .branch { meta, bam, bai, found ->
            existing: found != null
                return [ meta, found ]
            missing : true
                return [ meta, bam, bai ]
        }

    ch_counts_existing = ch_bam_branched.existing
        .map { meta, counts ->
            log.info "  [REUSE]    ${meta.id}  ->  ${counts.name}"
            return [ meta, counts ]
        }

    ch_bam_missing = ch_bam_branched.missing
        .map { meta, bam, bai ->
            log.info "  [GENERATE] ${meta.id}  ->  no precomputed counts, running CollectReadCounts"
            return [ meta, bam, bai ]
        }
        .combine(ch_intervals)
        .map { meta, bam, bai, intervals -> [ meta, bam, bai, intervals ] }

    GATK4_COLLECTREADCOUNTS (
        ch_bam_missing,
        ch_fasta,
        ch_fai,
        ch_dict
    )
    ch_versions = ch_versions.mix(GATK4_COLLECTREADCOUNTS.out.versions.ifEmpty([]))

    ch_counts_new = GATK4_COLLECTREADCOUNTS.out.tsv
        .mix(GATK4_COLLECTREADCOUNTS.out.hdf5)
        .filter { meta, f -> f != null }

    ch_counts = ch_counts_existing.mix(ch_counts_new)

    //
    // Aggregate all counts for cohort steps
    //
    ch_all_counts_files = ch_counts.map { meta, counts -> counts }.collect()
    ch_all_counts = ch_all_counts_files
        .map { counts -> [ [ id: 'cohort' ], counts ] }

    //
    // AnnotateIntervals (needed by FilterIntervals and CNVCaller)
    //
    // Optional tracks: mappability and segmental-duplication (recommended by GATK).
    // Each track may be plain .bed or bgzipped .bed.gz (then a .tbi must sit next to it).
    //
    ch_intervals_meta = ch_intervals.map { intervals -> [ [ id: 'intervals' ], intervals ] }
    ch_dummy = Channel.value( [ [ id: 'null' ], [] ] )

    def build_track_channel = { String track_param, String label ->
        if (!track_param) return [ ch_dummy, ch_dummy ]
        def f = file(track_param, checkIfExists: true)
        def tbi_file = file("${track_param}.tbi", checkIfExists: false)
        def tbi = tbi_file.exists() ? tbi_file : []
        log.info "AnnotateIntervals: using ${label} track ${f.name}" +
                 (tbi instanceof List ? " (no .tbi)" : " (+${tbi.name})")
        return [
            Channel.value( [ [ id: label ], f   ] ),
            Channel.value( [ [ id: label ], tbi ] )
        ]
    }

    def (ch_mappa, ch_mappa_tbi)   = build_track_channel(params.mappability_track,           'mappability')
    def (ch_segdup, ch_segdup_tbi) = build_track_channel(params.segmental_duplication_track, 'segdup')

    GATK4_ANNOTATEINTERVALS (
        ch_intervals_meta,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_mappa,
        ch_mappa_tbi,
        ch_segdup,
        ch_segdup_tbi
    )

    ch_annotated = GATK4_ANNOTATEINTERVALS.out.annotated_intervals
        .map { meta, annotated -> annotated }

    //
    // FilterIntervals (optional)
    //
    if (params.use_filter_intervals) {
        log.info "FilterIntervals ENABLED: noisy bins will be removed from PON intervals."

        ch_filter_input = GATK4_BEDTOINTERVALLIST.out.interval_list.first()
        ch_annotated_meta = GATK4_ANNOTATEINTERVALS.out.annotated_intervals.first()

        GATK4_FILTERINTERVALS (
            ch_filter_input,
            ch_annotated_meta,
            ch_all_counts_files
        )
        ch_versions = ch_versions.mix(GATK4_FILTERINTERVALS.out.versions.first())

        ch_filtered_intervals = GATK4_FILTERINTERVALS.out.filtered_intervals
            .map { meta, intervals -> intervals }
    } else {
        log.info "FilterIntervals DISABLED (--use_filter_intervals false): using unfiltered intervals."
        ch_filtered_intervals = ch_intervals
    }

    //
    // DetermineGermlineContigPloidy (COHORT)
    //
    ch_ploidy_input = ch_all_counts
        .combine(ch_filtered_intervals)
        .map { meta, counts, intervals -> [ meta, counts, intervals, [] ] }

    ch_ploidy_model_empty = Channel.value( [ [ id: 'ploidy_model' ], [] ] )

    GATK4_DETERMINEGERMLINECONTIGPLOIDY (
        ch_ploidy_input,
        ch_ploidy_model_empty,
        ch_ploidy_priors
    )
    ch_versions = ch_versions.mix(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.versions.first())

    ch_ploidy_calls = GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls
        .map { meta, calls -> calls }

    //
    // GermlineCNVCaller (COHORT)
    //
    ch_cnvcaller_input = ch_all_counts
        .combine(ch_filtered_intervals)
        .combine(ch_ploidy_calls)
        .combine(ch_annotated)
        .map { meta, counts, intervals, ploidy, annotated -> [ meta, counts, intervals, ploidy, [], annotated ] }

    GATK4_GERMLINECNVCALLER ( ch_cnvcaller_input )
    ch_versions = ch_versions.mix(GATK4_GERMLINECNVCALLER.out.versions.first())

    //
    // Build PON manifest
    //
    // Build per-sample (id, source) channel
    ch_sample_sources = ch_counts_existing.map { meta, c -> [ meta.id, 'reused'   ] }
        .mix(ch_counts_new.map     { meta, c -> [ meta.id, 'computed' ] })
        .toList()
        .map { entries ->
            def ids     = entries.collect { it[0] }
            def sources = entries.collectEntries { [ (it[0]): it[1] ] }
            return [ ids, sources ]
        }

    def filter_params_map = [
        minimum_gc_content        : params.minimum_gc_content,
        maximum_gc_content        : params.maximum_gc_content,
        minimum_mappability       : params.minimum_mappability,
        maximum_segmental_dup     : params.maximum_segmental_dup,
        low_count_filter_count    : params.low_count_filter_count_threshold,
        low_count_filter_pct      : params.low_count_filter_percentage_of_samples,
        extreme_count_min_pct     : params.extreme_count_filter_minimum_percentile,
        extreme_count_max_pct     : params.extreme_count_filter_maximum_percentile,
        extreme_count_pct_samples : params.extreme_count_filter_percentage_of_samples
    ]

    ch_bed_path   = ch_bed.map   { meta, b -> b }
    ch_fasta_path = ch_fasta.map { meta, f -> f }

    ch_manifest_in = ch_sample_sources
        .combine(ch_bed_path)
        .combine(ch_fasta_path)
        .combine(ch_intervals)
        .combine(ch_filtered_intervals)

    PON_MANIFEST (
        ch_manifest_in.map { ids, sources, bed, fasta, iv, fi -> ids },
        ch_manifest_in.map { ids, sources, bed, fasta, iv, fi -> sources },
        ch_manifest_in.map { ids, sources, bed, fasta, iv, fi -> bed },
        ch_manifest_in.map { ids, sources, bed, fasta, iv, fi -> fasta },
        ch_manifest_in.map { ids, sources, bed, fasta, iv, fi -> iv },
        ch_manifest_in.map { ids, sources, bed, fasta, iv, fi -> fi },
        params.use_filter_intervals,
        workflow.manifest.version,
        filter_params_map,
        params.bed_padding ?: 0
    )
    ch_versions = ch_versions.mix(PON_MANIFEST.out.versions)

    emit:
    model              = GATK4_GERMLINECNVCALLER.out.cohortmodel
    cohort_calls       = GATK4_GERMLINECNVCALLER.out.cohortcalls
    ploidy_calls       = GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls
    ploidy_model       = GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.model
    intervals          = ch_intervals
    filtered_intervals = ch_filtered_intervals
    counts             = ch_counts
    annotated          = ch_annotated
    manifest           = PON_MANIFEST.out.manifest
    versions           = ch_versions
}
