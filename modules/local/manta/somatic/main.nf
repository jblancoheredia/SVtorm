process MANTA_SOMATIC {
    tag "$meta.patient"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/manta:1.6.0--h9ee0642_1' :
        'quay.io/biocontainers/manta:1.6.0--h9ee0642_1' }"

    input:
    tuple val(meta) , path(tumour_bam), path(tumour_bai), path(normal_bam), path(normal_bai)
    tuple val(meta2), path(target_bed) 
    tuple val(meta3), path(bed_index)
    tuple val(meta4), path(fasta)
    tuple val(meta5), path(fai)
    path(config)

    output:
    tuple val(meta), path("*.candidate_small_indels.vcf.gz")     , emit: candidate_small_indels_vcf
    tuple val(meta), path("*.candidate_small_indels.vcf.gz.tbi") , emit: candidate_small_indels_vcf_tbi
    tuple val(meta), path("*.candidate_sv.vcf.gz")               , emit: candidate_sv_vcf
    tuple val(meta), path("*.candidate_sv.vcf.gz.tbi")           , emit: candidate_sv_vcf_tbi
    tuple val(meta), path("*.diploid_sv.vcf.gz")                 , emit: diploid_sv_vcf
    tuple val(meta), path("*.diploid_sv.vcf.gz.tbi")             , emit: diploid_sv_vcf_tbi
    tuple val(meta), path("*.somatic_sv.vcf.gz")                 , emit: somatic_sv_vcf
    tuple val(meta), path("*.somatic_sv.vcf.gz.tbi")             , emit: somatic_sv_vcf_tbi
    tuple val(meta), path("*.manta.unfiltered.vcf")              , emit: vcf
    tuple val(meta), path("*.tsv")                               , emit: metrics
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    def options_manta = target_bed ? "--callRegions $target_bed" : ""
    def config_option = config ? "--config ${config}" : ""
    """
    configManta.py \\
        --tumorBam $tumour_bam \\
        --normalBam $normal_bam \\
        --reference $fasta \\
        ${config_option} \\
        --runDir manta \\
        $options_manta \\
        $args

    python manta/runWorkflow.py -m local -j $task.cpus

    mv manta/results/variants/candidateSmallIndels.vcf.gz \\
        ${prefix}.candidate_small_indels.vcf.gz
    mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \\
        ${prefix}.candidate_small_indels.vcf.gz.tbi
    mv manta/results/variants/candidateSV.vcf.gz \\
        ${prefix}.candidate_sv.vcf.gz
    mv manta/results/variants/candidateSV.vcf.gz.tbi \\
        ${prefix}.candidate_sv.vcf.gz.tbi
    mv manta/results/variants/diploidSV.vcf.gz \\
        ${prefix}.diploid_sv.vcf.gz
    mv manta/results/variants/diploidSV.vcf.gz.tbi \\
        ${prefix}.diploid_sv.vcf.gz.tbi
    mv manta/results/variants/somaticSV.vcf.gz \\
        ${prefix}.somatic_sv.vcf.gz
    mv manta/results/variants/somaticSV.vcf.gz.tbi \\
        ${prefix}.somatic_sv.vcf.gz.tbi
    mv manta/results/stats/alignmentStatsSummary.txt \\
        ${prefix}.alignmentStatsSummary.tsv
    mv manta/results/stats/svCandidateGenerationStats.tsv \\
        ${prefix}.svCandidateGenerationStats.tsv
    mv manta/results/stats/svLocusGraphStats.tsv \\
        ${prefix}.svLocusGraphStats.tsv

    cp ${prefix}.candidate_sv.vcf.gz ${prefix}.manta.unfiltered.vcf.gz

    gunzip ${prefix}.manta.unfiltered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$( configManta.py --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    echo "" | gzip > ${prefix}.candidate_small_indels.vcf.gz
    touch ${prefix}.candidate_small_indels.vcf.gz.tbi
    echo "" | gzip > ${prefix}.candidate_sv.vcf.gz
    touch ${prefix}.candidate_sv.vcf.gz.tbi
    echo "" | gzip > ${prefix}.diploid_sv.vcf.gz
    touch ${prefix}.diploid_sv.vcf.gz.tbi
    echo "" | gzip > ${prefix}.somatic_sv.vcf.gz
    touch ${prefix}.somatic_sv.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$( configManta.py --version )
    END_VERSIONS
    """
}
