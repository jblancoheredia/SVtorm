process IGVREPORTS {
    tag "$meta.patient"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/igv-reports:1.12.0--pyh7cba7a3_0':
        'biocontainers/igv-reports:1.12.0--pyh7cba7a3_0' }"

    input:
    tuple val(meta) , path(tumour_bam), path(tumour_bai), path(normal_bam), path(normal_bai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta3), path(sites)

    output:
    tuple val(meta), path("*.html"), emit: report
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    def fasta = fasta ? "--fasta ${fasta}" : ""
    def track_arg = tracks ? "--tracks ${tumour_bam}" : ""
    """
    create_report \\
        ${sites} \\
        --fasta ${fasta} \\
        --tracks ${tumour_bam} \\
        --output ${prefix}_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igvreports: \$(python -c "import igv_reports; print(igv_reports.__version__)")
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    touch ${prefix}_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igvreports: \$(python -c "import igv_reports; print(igv_reports.__version__)")
    END_VERSIONS
    """
}
