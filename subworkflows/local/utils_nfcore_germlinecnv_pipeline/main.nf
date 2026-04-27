/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UTILS_NFCORE_GERMLINECNV_PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Utility subworkflows for pipeline initialization and completion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PIPELINE_INITIALISATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args // array: List of positional nextflow CLI args
    outdir            // string: The output directory where the results will be saved

    main:

    //
    // Print version and exit if required
    //
    if (version) {
        log.info "${workflow.manifest.name} ${workflow.manifest.version}"
        System.exit(0)
    }

    //
    // Print help message if required
    //
    if (help) {
        log.info helpMessage()
        System.exit(0)
    }

    //
    // Print parameter summary
    //
    log.info paramsSummary()

    //
    // Create samplesheet channel
    //
    ch_samplesheet = Channel.empty()
    
    if (params.input) {
        ch_samplesheet = Channel
            .fromPath(params.input)
            .splitCsv(header: true, sep: ',')
            .map { row ->
                def meta = [:]
                meta.id = row.sample
                meta.single_end = false
                
                def bam = file(row.bam, checkIfExists: true)
                def bai = row.bai ? file(row.bai, checkIfExists: true) : file("${row.bam}.bai", checkIfExists: false)
                if (!bai.exists()) {
                    bai = file(row.bam.replaceAll(/\.bam$/, '.bai'), checkIfExists: false)
                }
                
                return [ meta, bam, bai ]
            }
            .view { "Parsed sample: ${it[0].id}" }
    }

    emit:
    samplesheet = ch_samplesheet
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PIPELINE_COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email             // string: email address
    email_on_fail     // string: email address on failure
    plaintext_email   // boolean: Use plaintext email
    outdir            // string: Output directory
    monochrome_logs   // boolean: Disable coloured log outputs
    hook_url          // string: Hook URL for notifications

    main:

    //
    // Completion summary
    //
    workflow.onComplete {
        completionSummary(monochrome_logs)
        
        if (email || email_on_fail) {
            completionEmail(email, email_on_fail, plaintext_email, outdir)
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Help message
//
def helpMessage() {
    return """
    =========================================
    nf-core/germlinecnv v${workflow.manifest.version}
    =========================================
    
    Usage:
        nextflow run main.nf --input samplesheet.csv --fasta ref.fa --bed targets.bed --outdir results -profile docker

    Mandatory arguments:
        --input         Path to samplesheet CSV
        --fasta         Path to reference FASTA
        --bed           Path to BED file with target regions
        --outdir        Output directory

    Pipeline modes:
        --mode          Execution mode: 'pon', 'case', or 'full' (default: full)

    For CASE mode:
        --pon_model     Path to existing PON model
        --ploidy_model  Path to existing ploidy model

    Optional:
        --mane_file     Path to MANE annotation file
        -profile        Configuration profile (docker, singularity, conda)
    """.stripIndent()
}

//
// Parameter summary
//
def paramsSummary() {
    return """
    =========================================
    nf-core/germlinecnv v${workflow.manifest.version}
    =========================================
    input       : ${params.input}
    mode        : ${params.mode}
    fasta       : ${params.fasta}
    bed         : ${params.bed}
    outdir      : ${params.outdir}
    -----------------------------------------
    """.stripIndent()
}

//
// Completion email
//
def completionEmail(email, email_on_fail, plaintext_email, outdir) {
    def subject = "[nf-core/germlinecnv] Successful: ${workflow.runName}"
    if (!workflow.success) {
        subject = "[nf-core/germlinecnv] FAILED: ${workflow.runName}"
    }

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at : ${workflow.complete}
        Duration     : ${workflow.duration}
        Success      : ${workflow.success}
        workDir      : ${workflow.workDir}
        exit status  : ${workflow.exitStatus}
        """.stripIndent()

    if (email && workflow.success) {
        sendMail(to: email, subject: subject, body: msg)
    }
    if (email_on_fail && !workflow.success) {
        sendMail(to: email_on_fail, subject: subject, body: msg)
    }
}

//
// Completion summary
//
def completionSummary(monochrome_logs) {
    Map colors = logColours(monochrome_logs)
    if (workflow.success) {
        log.info "-${colors.purple}[nf-core/germlinecnv]${colors.green} Pipeline completed successfully${colors.reset}-"
    } else {
        log.info "-${colors.purple}[nf-core/germlinecnv]${colors.red} Pipeline completed with errors${colors.reset}-"
    }
}

//
// Log colours
//
def logColours(monochrome_logs) {
    Map colorcodes = [:]
    colorcodes['reset']  = monochrome_logs ? '' : "\033[0m"
    colorcodes['red']    = monochrome_logs ? '' : "\033[0;31m"
    colorcodes['green']  = monochrome_logs ? '' : "\033[0;32m"
    colorcodes['yellow'] = monochrome_logs ? '' : "\033[0;33m"
    colorcodes['blue']   = monochrome_logs ? '' : "\033[0;34m"
    colorcodes['purple'] = monochrome_logs ? '' : "\033[0;35m"
    colorcodes['cyan']   = monochrome_logs ? '' : "\033[0;36m"
    colorcodes['white']  = monochrome_logs ? '' : "\033[0;37m"
    return colorcodes
}
