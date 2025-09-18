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
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fasta_fai)
    tuple val(meta5), path(known_sites)
    tuple val(meta6), path(known_sites_tbi)
    path(refflat)
    path(bed)
    path(blocklist)
    path(bwa_index)
    path(kraken2db)
    path(pon_dir)

    output:
    tuple val(meta), path("*.recall.all_calls_avk.vcf"), emit: vcf
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def task_cpus = task.cpus <= 8 ? task.cpus : 8
    def prefix = task.ext.prefix ?: "${meta.patient}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
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

    CollectGridssMetrics \\
        INPUT=${prefix}_N_filtered.bam \\
        THRESHOLD_COVERAGE=100000 \\
        PROGRAM=MeanQualityByCycle \\
        PROGRAM=CollectGcBiasMetrics \\
        PROGRAM=CollectInsertSizeMetrics \\
        PROGRAM=QualityScoreDistribution  \\
        PROGRAM=CollectQualityYieldMetrics \\
        PROGRAM=CollectBaseDistributionByCycle \\
        PROGRAM=CollectAlignmentSummaryMetrics  \\
        PROGRAM=CollectSequencingArtifactMetrics \\
        OUTPUT=${prefix}-N.bam.gridss.working/${prefix}-N.bam \\
        INTERVALS=${interval_list} \\
        REF_FLAT=${refflat} \\
        REFERENCE_SEQUENCE=${fasta} \\
        DB_SNP=${known_sites}

    gridss_extract_overlapping_fragments \\
        --targetbed ${prefix}.interval_list.bed \\
        -t ${task_cpus} \\
        -o ${prefix}-N.bam \\
        ${prefix}_N_filtered.bam

    mkdir ${prefix}-T.bam.gridss.working

    CollectGridssMetrics \\
        INPUT=${prefix}_T_filtered.bam \\
        THRESHOLD_COVERAGE=100000 \\
        PROGRAM=MeanQualityByCycle \\
        PROGRAM=CollectGcBiasMetrics \\
        PROGRAM=CollectInsertSizeMetrics \\
        PROGRAM=QualityScoreDistribution  \\
        PROGRAM=CollectQualityYieldMetrics \\
        PROGRAM=CollectBaseDistributionByCycle \\
        PROGRAM=CollectAlignmentSummaryMetrics  \\
        PROGRAM=CollectSequencingArtifactMetrics \\
        OUTPUT=${prefix}-T.bam.gridss.working/${prefix}-T.bam \\
        INTERVALS=${interval_list} \\
        VALIDATION_STRINGENCY=LENIENT \\
        REF_FLAT=${refflat} \\
        REFERENCE_SEQUENCE=${fasta} \\
        DB_SNP=${known_sites}

    gridss_extract_overlapping_fragments \\
        --targetbed ${prefix}.interval_list.bed \\
        -t ${task_cpus} \\
        -o ${prefix}-T.bam \\
        ${prefix}_T_filtered.bam

    gridss \\
        --threads ${task_cpus} \\
        --labels "NORMAL",${prefix} \\
        --jvmheap ${task.memory.toGiga() - 1}g \\
        --otherjvmheap ${task.memory.toGiga() - 1}g \\
        --output ${prefix}_all_calls.vcf \\
        --reference ${fasta} \\
        --picardoptions VALIDATION_STRINGENCY=LENIENT \\
        -b ${blocklist} \\
        ${prefix}-N.bam \\
        ${prefix}-T.bam

    gridss_annotate_vcf_kraken2 \\
        -t ${task_cpus} \\
        -o ${prefix}.recall.all_calls_avk.vcf \\
        --kraken2db ${kraken2db} \\
        ${prefix}_all_calls.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.recall.all_calls_avk.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
}
