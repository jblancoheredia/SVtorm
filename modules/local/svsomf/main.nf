process SVSOMF {
    tag "$meta.patient"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/svsomf:1.4.0':
        'blancojmskcc/svsomf:1.4.0' }"

    input:
    tuple val(meta) , path(vcf)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    path(pon_dir)

    output:
    tuple val(meta), path("*.unfiltered.vcf"), emit: vcf
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    def bam = "${tumour_bam}"
    """
    SVsomF \\
        --pondir ${pon_dir} \\
        --ref ${fasta} \\
        --input ${vcf} \\
        --output ${prefix}.recall.filtered.vcf \\
        --fulloutput ${prefix}.recall.unfiltered.vcf \\
        --plotdir . \\
        --normalordinal 1 \\
        --tumourordinal 2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svsomf: "1.4.0"
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    touch ${prefix}.recall.unfiltered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drawsv: "6.0.6"
    END_VERSIONS
    """
}
