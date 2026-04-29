process VALIDATE_PON_MANIFEST {
    tag "validate_pon"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/coreutils:9.1'
        : 'quay.io/biocontainers/coreutils:9.1'}"

    input:
    path manifest    // pon_manifest.json from a previous PON run
    path bed         // current case BED
    path fasta       // current case FASTA

    output:
    path "validation_report.txt", emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def strict = (params.strict_pon_validation == null) ? true : params.strict_pon_validation
    """
    set -euo pipefail

    bed_md5=\$(md5sum ${bed}   | awk '{print \$1}')
    fa_md5=\$(md5sum ${fasta}  | awk '{print \$1}')

    pon_bed_md5=\$(grep -oE '"bed_md5"[[:space:]]*:[[:space:]]*"[^"]+"' ${manifest}   | head -1 | sed 's/.*"\\([^"]*\\)"\$/\\1/')
    pon_fa_md5=\$(grep -oE '"fasta_md5"[[:space:]]*:[[:space:]]*"[^"]+"' ${manifest} | head -1 | sed 's/.*"\\([^"]*\\)"\$/\\1/')

    {
        echo "=== PON manifest validation ==="
        echo "PON BED md5   : \${pon_bed_md5}"
        echo "Case BED md5  : \${bed_md5}"
        echo "PON FASTA md5 : \${pon_fa_md5}"
        echo "Case FASTA md5: \${fa_md5}"
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
