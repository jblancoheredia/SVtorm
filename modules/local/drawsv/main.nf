process DRAWSV {
    tag "$meta.patient"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/drawsv:1.0.2':
        'blancojmskcc/drawsv:1.0.2' }"

    input:
    tuple val(meta) , path(tumour_bam), path(tumour_bai), path(normal_bam), path(normal_bai)
    tuple val(meta2), path(tsv)
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
    def prefix = task.ext.prefix ?: "${meta.patient}"
    def bam = "${tumour_bam}"
    """
    cp "${workflow.projectDir}/bin/PRE_DRAWSVs.py" .
    cp "${workflow.projectDir}/bin/DrawSVs.R" .

    #python3 PRE_DRAWSVs.py \\
    PreDrawSVs \\
        --sample ${prefix} \\
        --input ${tsv} \\
        --genome ${genome} \\
        --annotations ${gtf} \\
        ${args}

    #Rscript DrawSVs.R \\
    Rscript /usr/local/bin/DrawSVs.R \\
        --SVs=${prefix}_DrawSVs.tsv \\
        --alignments=${bam} \\
        --annotation=${gtf}   \\
        --cytobands=${cytobands} \\
        --output=${prefix}_DrawSVs.pdf \\
        --transcriptSelection=canonical \\
        --minConfidenceForCircosPlot=High \\
        --proteinDomains=${protein_domains} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drawsv: "1.0.1"
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    touch ${prefix}_DrawSVs.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drawsv: "1.0.1"
    END_VERSIONS
    """
}
