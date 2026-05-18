process PON_MANIFEST {
    tag "pon_manifest"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:22.04'
        : 'docker.io/ubuntu:22.04'}"

    input:
    val  samples_used        // list of sample IDs
    val  counts_sources      // map: id -> "computed" | "reused"
    path bed                 // capture BED
    path fasta               // reference FASTA
    path intervals           // preprocessed interval_list
    path filtered_intervals  // filtered interval_list (or NO_FILE)
    val  filter_enabled
    val  pipeline_version
    val  filter_params       // map of filter thresholds (for documentation)
    val  bed_padding         // BedToIntervalList --PADDING value (0 if disabled)

    output:
    path "pon_manifest.json", emit: manifest
    path "pon_samples.tsv",   emit: samples_tsv
    path "versions.yml",      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def created_at = new java.util.Date().format("yyyy-MM-dd'T'HH:mm:ssXXX")
    def samples_json = groovy.json.JsonOutput.toJson(
        samples_used.collect { id ->
            [ id: id, counts_source: (counts_sources[id] ?: 'unknown') ]
        }
    )
    def filter_json   = groovy.json.JsonOutput.toJson(filter_params ?: [:])
    def filtered_name = filtered_intervals ? filtered_intervals.name : ''
    def tsv_rows = samples_used.collect { id ->
        "${id}\\t${counts_sources[id] ?: 'unknown'}"
    }.join('\\n')
    """
    set -euo pipefail

    bed_md5=\$(md5sum ${bed}      | awk '{print \$1}')
    fa_md5=\$(md5sum ${fasta}     | awk '{print \$1}')
    iv_md5=\$(md5sum ${intervals} | awk '{print \$1}')
    fi_md5=""
    if [ -n "${filtered_name}" ] && [ -s "${filtered_name}" ]; then
        fi_md5=\$(md5sum ${filtered_name} | awk '{print \$1}')
    fi

    cat > pon_manifest.json <<EOF
{
  "pipeline_version": "${pipeline_version}",
  "created_at": "${created_at}",
  "reference": {
    "fasta_path": "${fasta}",
    "fasta_md5": "\${fa_md5}"
  },
  "intervals": {
    "bed_path": "${bed}",
    "bed_md5": "\${bed_md5}",
    "bed_padding": ${bed_padding ?: 0},
    "preprocessed_intervals": "${intervals}",
    "preprocessed_intervals_md5": "\${iv_md5}",
    "filtered_intervals": "${filtered_name}",
    "filtered_intervals_md5": "\${fi_md5}"
  },
  "filter_intervals": {
    "enabled": ${filter_enabled ? 'true' : 'false'},
    "params": ${filter_json}
  },
  "samples": ${samples_json}
}
EOF

    # Sample TSV (id, counts_source)
    {
      printf "sample_id\\tcounts_source\\n"
      printf "%s\\n" "${tsv_rows}"
    } > pon_samples.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(md5sum --version | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch pon_manifest.json pon_samples.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: stub
    END_VERSIONS
    """
}
