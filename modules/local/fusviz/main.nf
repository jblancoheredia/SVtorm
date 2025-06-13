process FUSVIZ {
    tag "$meta.patient_id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/fusviz:1.0.2':
        'blancojmskcc/fusviz:1.0.2' }"

    input:
    tuple val(patient_id), 
          val(meta_normal), path(normal_bam), path(normal_bai),
          val(meta_tumour), path(tumour_bam), path(tumour_bai)
    tuple val(meta), path(tsv)
    path(gtf)
    val(genome)
    path(cytobands)
    path(protein_domains)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient_id}"
    def bam = "${tumour_bam}"
    """
    preFusViz \\
        --sample ${prefix} \\
        --input ${tsv} \\
        --genome ${genome} \\
        --annotations ${gtf} \\
        ${args}

    FusViz \\
        --SVs=${prefix}_FusViz.tsv \\
        --alignments=${bam} \\
        --annotation=${gtf}   \\
        --cytobands=${cytobands} \\
        --output=${prefix}_FusViz.pdf \\
        --transcriptSelection=canonical \\
        --minConfidenceForCircosPlot=High \\
        --proteinDomains=${protein_domains} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusviz: "1.0.2"
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient_id}"
    """
    touch ${prefix}_FusViz.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusviz: "1.0.2"
    END_VERSIONS
    """
}
