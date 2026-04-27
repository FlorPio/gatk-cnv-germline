process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0' }"

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env python

    import yaml
    import platform
    from textwrap import dedent

    def _make_versions_html(versions):
        html = [
            dedent(
                '''\\
                <style>
                #nf-core-versions tbody:nth-child(even) {
                    background-color: #f2f2f2;
                }
                </style>
                <table class="table" style="width:100%" id="nf-core-versions">
                    <thead>
                        <tr>
                            <th> Process Name </th>
                            <th> Software </th>
                            <th> Version </th>
                        </tr>
                    </thead>
                '''
            )
        ]
        for process, tmp_versions in sorted(versions.items()):
            html.append("<tbody>")
            for i, (tool, version) in enumerate(sorted(tmp_versions.items())):
                html.append(
                    dedent(
                        f'''\\
                        <tr>
                            <td><samp>{process if (i == 0) else ''}</samp></td>
                            <td><samp>{tool}</samp></td>
                            <td><samp>{version}</samp></td>
                        </tr>
                        '''
                    )
                )
            html.append("</tbody>")
        html.append("</table>")
        return "\\n".join(html)

    # Load versions
    with open("${versions}") as f:
        versions = yaml.safe_load(f) or {}

    # Add workflow versions
    versions["Workflow"] = {
        "Nextflow": "$workflow.nextflow.version",
        "${workflow.manifest.name}": "${workflow.manifest.version}"
    }

    # Write plain YAML
    with open("software_versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    # Write MultiQC YAML
    mqc_yml = {
        "id": "software_versions",
        "section_name": "${workflow.manifest.name} Software Versions",
        "section_href": "https://github.com/nf-core/germlinecnv",
        "plot_type": "html",
        "description": "Software versions collected at run time from the software output.",
        "data": _make_versions_html(versions)
    }
    
    with open("software_versions_mqc.yml", "w") as f:
        yaml.dump(mqc_yml, f, default_flow_style=False)

    # Write versions.yml
    versions_yml = {
        "${task.process}": {
            "python": platform.python_version(),
            "yaml": yaml.__version__
        }
    }
    
    with open("versions.yml", "w") as f:
        yaml.dump(versions_yml, f, default_flow_style=False)
    """

    stub:
    """
    touch software_versions.yml
    touch software_versions_mqc.yml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
