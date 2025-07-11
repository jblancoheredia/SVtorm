process SURVIVOR_STATS {
    tag "$meta.patient"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/survivor:1.0.7--h9a82719_1':
        'biocontainers/survivor:1.0.7--h9a82719_1' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val(minsv)          // Min SV size (-1 to disable)
    val(maxsv)          // Max SV size (-1 to disable)
    val(minnumreads)    // Min number of reads support: RE flag (-1 to disable)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    def is_compressed = vcf.getName().endsWith(".gz") ? true : false
    vcf_name = vcf.getName().replace(".gz", "")

    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $vcf > $vcf_name
    fi

    SURVIVOR \\
        stats \\
        $vcf_name \\
        $minsv \\
        $maxsv \\
        $minnumreads \\
        ${prefix}.stats

    rm $vcf_name

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.patient}"

    """
    touch ${prefix}.stats
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """
}
