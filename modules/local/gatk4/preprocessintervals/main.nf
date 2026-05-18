process GATK4_PREPROCESSINTERVALS {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::gatk4=4.5.0.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta),  path(intervals)        // BED or interval_list
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(dict)
    tuple val(meta5), path(exclude_intervals) // optional, can be []

    output:
    tuple val(meta), path("*.interval_list"), emit: interval_list
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def exclude_arg = (exclude_intervals && !(exclude_intervals instanceof List && exclude_intervals.isEmpty()))
                        ? "--exclude-intervals ${exclude_intervals}"
                        : ''
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        PreprocessIntervals \\
        --intervals ${intervals} \\
        --reference ${fasta} \\
        --sequence-dictionary ${dict} \\
        ${exclude_arg} \\
        --output ${prefix}.interval_list \\
        --tmp-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/.*GATK *v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.interval_list
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: stub
    END_VERSIONS
    """
}
