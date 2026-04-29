process GATK4_FILTERINTERVALS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(intervals)
    tuple val(meta2), path(annotated_intervals)
    path(read_counts)

    output:
    tuple val(meta), path("*.interval_list"), emit: filtered_intervals
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_filtered"
    def annotated_arg = annotated_intervals ? "--annotated-intervals ${annotated_intervals}" : ''
    def counts_args   = read_counts ? read_counts.collect { "--input ${it}" }.join(' ') : ''

    def avail_mem = 3072
    if (task.memory) {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        FilterIntervals \\
        --intervals ${intervals} \\
        ${annotated_arg} \\
        ${counts_args} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --output ${prefix}.interval_list \\
        --tmp-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n '/GATK.*v/s/.*v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_filtered"
    """
    touch ${prefix}.interval_list
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n '/GATK.*v/s/.*v//p')
    END_VERSIONS
    """
}
