process VALIDATE_PON_MANIFEST {
    tag "validate_pon"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:22.04'
        : 'docker.io/ubuntu:22.04'}"

    input:
    path manifest    // pon_manifest.json from a previous PON run
    path bed         // current case BED
    path fasta       // current case FASTA

    output:
    path "validation_report.txt", emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def strict       = (params.strict_pon_validation == null) ? true : params.strict_pon_validation
    def case_padding = params.bed_padding ?: 0
    """
    set -euo pipefail

    bed_md5=\$(md5sum ${bed}   | awk '{print \$1}')
    fa_md5=\$(md5sum ${fasta}  | awk '{print \$1}')

    pon_bed_md5=\$(grep -oE '"bed_md5"[[:space:]]*:[[:space:]]*"[^"]+"' ${manifest}   | head -1 | sed 's/.*"\\([^"]*\\)"\$/\\1/')
    pon_fa_md5=\$(grep -oE '"fasta_md5"[[:space:]]*:[[:space:]]*"[^"]+"' ${manifest} | head -1 | sed 's/.*"\\([^"]*\\)"\$/\\1/')
    pon_padding=\$(grep -oE '"bed_padding"[[:space:]]*:[[:space:]]*[0-9]+' ${manifest} | head -1 | grep -oE '[0-9]+\$' || echo "0")

    {
        echo "=== PON manifest validation ==="
        echo "PON BED md5    : \${pon_bed_md5}"
        echo "Case BED md5   : \${bed_md5}"
        echo "PON FASTA md5  : \${pon_fa_md5}"
        echo "Case FASTA md5 : \${fa_md5}"
        echo "PON bed_padding: \${pon_padding}"
        echo "Case bed_padding: ${case_padding}"
    } > validation_report.txt

    fail=0
    if [ "\${pon_bed_md5}" != "\${bed_md5}" ]; then
        echo "MISMATCH: BED differs between PON and case." >> validation_report.txt
        fail=1
    fi
    if [ "\${pon_fa_md5}" != "\${fa_md5}" ]; then
        echo "MISMATCH: FASTA differs between PON and case." >> validation_report.txt
        fail=1
    fi
    if [ "\${pon_padding}" != "${case_padding}" ]; then
        echo "MISMATCH: bed_padding differs (PON=\${pon_padding}, case=${case_padding}). The intervals will not align." >> validation_report.txt
        fail=1
    fi

    if [ "\$fail" -eq 1 ]; then
        if [ "${strict}" = "true" ]; then
            cat validation_report.txt 1>&2
            echo "ERROR: PON/case BED or FASTA mismatch. Set --strict_pon_validation false to override." 1>&2
            exit 1
        else
            echo "WARNING: differences detected but --strict_pon_validation=false; continuing." >> validation_report.txt
        fi
    else
        echo "OK: PON and case share BED and FASTA." >> validation_report.txt
    fi
    """

    stub:
    """
    echo "stub" > validation_report.txt
    """
}
