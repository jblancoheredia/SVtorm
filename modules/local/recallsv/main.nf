process RECALL_SV {
    tag "$meta.patient"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/gridss:2.13.2':
        'blancojmskcc/gridss:2.13.2' }"

    input:
    tuple val(meta), 
          val(meta0), path(tumour_bam), path(tumour_bai),
          val(meta1), path(normal_bam), path(normal_bai),
          val(meta2), path(interval_list)
    tuple val(meta3), path(knownsites_i)
    tuple val(meta4), path(known_sites)
    tuple val(meta5), path(fasta_fai)
    tuple val(meta6), path(fasta)
    path(blocklist)
    path(bwa_index)
    path(kraken2db)
    path(pon_dir)
    path(refflat)
    path(bed)

    output:
    tuple val(meta), path("*.recall.all_calls_avk.vcf"), emit: vcf
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    def VERSION = '2.13.2'
    def bwa = bwa_index ? "cp -s ${bwa_index}/* ." : ""
    """
    source activate gridss

    samtools view -@ ${task_cpus} -h -F 256 -o ${prefix}_N_filtered.bam ${normal_bam}
    samtools index -@ ${task_cpus} ${prefix}_N_filtered.bam

    samtools view -@ ${task_cpus} -h -F 256 -o ${prefix}_T_filtered.bam ${tumour_bam}
    samtools index -@ ${task_cpus} ${prefix}_T_filtered.bam

    rm ${fasta} ${fasta_fai}

    ${bwa}

    grep -v "@" ${interval_list} > ${prefix}.interval_list.bed

    mkdir ${prefix}-N.bam.gridss.working

    CollectGridssMetrics  \\
        REF_FLAT=${refflat} \\
        DB_SNP=${known_sites} \\
        THRESHOLD_COVERAGE=100000 \\
        INTERVALS=${interval_list} \\
        PROGRAM=MeanQualityByCycle  \\
        REFERENCE_SEQUENCE=${fasta}  \\
        PROGRAM=CollectGcBiasMetrics  \\
        INPUT=${prefix}_N_filtered.bam \\
        PROGRAM=CollectInsertSizeMetrics \\
        PROGRAM=QualityScoreDistribution  \\
        PROGRAM=CollectQualityYieldMetrics \\
        PROGRAM=CollectBaseDistributionByCycle \\
        PROGRAM=CollectAlignmentSummaryMetrics  \\
        PROGRAM=CollectSequencingArtifactMetrics \\
        OUTPUT=${prefix}-N.bam.gridss.working/${prefix}-N.bam

    gridss_extract_overlapping_fragments \\
        -t 8 \\
        -o ${prefix}-N.bam \\
        --targetbed ${prefix}.interval_list.bed \\
        ${prefix}_N_filtered.bam

    mkdir ${prefix}-T.bam.gridss.working

    CollectGridssMetrics  \\
        REF_FLAT=${refflat} \\
        DB_SNP=${known_sites} \\
        THRESHOLD_COVERAGE=100000 \\
        INTERVALS=${interval_list} \\
        PROGRAM=MeanQualityByCycle  \\
        REFERENCE_SEQUENCE=${fasta}  \\
        PROGRAM=CollectGcBiasMetrics  \\
        INPUT=${prefix}_T_filtered.bam \\
        PROGRAM=CollectInsertSizeMetrics \\
        PROGRAM=QualityScoreDistribution  \\
        PROGRAM=CollectQualityYieldMetrics \\
        PROGRAM=CollectBaseDistributionByCycle \\
        PROGRAM=CollectAlignmentSummaryMetrics  \\
        PROGRAM=CollectSequencingArtifactMetrics \\
        OUTPUT=${prefix}-T.bam.gridss.working/${prefix}-N.bam

    gridss_extract_overlapping_fragments \\
        -t 8 \\
        -o ${prefix}-T.bam \\
        --targetbed ${prefix}.interval_list.bed \\
        ${prefix}_T_filtered.bam

    gridss \\
        --threads 8 \\
        -b ${blocklist} \\
        --reference ${fasta} \\
        --labels "NORMAL",${prefix} \\
        --output ${prefix}_all_calls.vcf \\
        --jvmheap ${task.memory.toGiga() - 1}g \\
        --otherjvmheap ${task.memory.toGiga() - 1}g \\
        --picardoptions VALIDATION_STRINGENCY=LENIENT  \\
        ${prefix}-N.bam \\
        ${prefix}-T.bam

    gridss_annotate_vcf_kraken2 \\
        -t 8 \\
        --kraken2db ${kraken2db} \\
        -o ${prefix}.recall.all_calls_avk.vcf \\
        ${prefix}_all_calls.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    def VERSION = '2.13.2'
    """
    touch ${prefix}.recall.all_calls_avk.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
}
