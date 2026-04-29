process ANNOTATE_CNV_EXONS {
    tag "${meta.id}"
    label 'process_single'

    container "docker.io/florpio/cnv-annotate-r:1.0"

    input:
    tuple val(meta), path(vcf)
    path(mane_annotation)
    path(bed)
    path(rscript)
    path(genes_list)

    output:
    tuple val(meta), path("*.txt"),      emit: annotated
    tuple val(meta), path("*.tsv"),      emit: annotated_tsv, optional: true
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genes_arg = genes_list ? "--genes ${genes_list}" : ''
    def bed_arg   = bed        ? "--bed ${bed}"          : ''
    """
    # Decompress VCF if gzipped
    if [[ "${vcf}" == *.gz ]]; then
        gunzip -c ${vcf} > ${prefix}_segments.vcf
        VCF_INPUT="${prefix}_segments.vcf"
    else
        VCF_INPUT="${vcf}"
    fi

    Rscript ${rscript} \\
        --vcf \${VCF_INPUT} \\
        --MANE ${mane_annotation} \\
        ${bed_arg} \\
        ${genes_arg} \\
        ${args}

    # Rename output to include sample prefix if needed
    for f in *.txt; do
        if [ -f "\$f" ] && [[ "\$f" != "${prefix}"* ]]; then
            mv "\$f" "${prefix}_annotated.txt" 2>/dev/null || true
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_annotated.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
    END_VERSIONS
    """
}
