process DRAWSV {
    tag "$meta.patient"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/svtorm_drawsv:7.1.3':
        'blancojmskcc/svtorm_drawsv:7.1.3' }"

    input:
    tuple val(meta), 
          val(meta0), path(tumour_bam), path(tumour_bai),
          val(meta1), path(normal_bam), path(normal_bai),
          val(meta2), path(tsv)
    path(gtf)
    path(cytobands)
    path(chromosomes)
    path(protein_domains)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    def bam = "${tumour_bam}"
    """
    touch ${prefix}_DrawSV.pdf
    
    DrawSV \\
        --SVs ${tsv} \\
        --alignments ${bam} \\
        --annotation ${gtf}   \\
        --cytobands ${cytobands} \\
        --chromosomes ${chromosomes} \\
        --output ${prefix}_DrawSV.pdf \\
        --proteinDomains=${protein_domains}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drawsv: "7.1.3"
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    touch ${prefix}_DrawSV.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drawsv: "7.1.3"
    END_VERSIONS
    """
}
