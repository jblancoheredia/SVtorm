process SVABA {
    tag "$meta.patient"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/svaba:1.2.0' :
        'blancojmskcc/svaba:1.2.0' }"

    input:
    tuple val(meta) , path(tumour_bam), path(tumour_bai), path(normal_bam), path(normal_bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    path(dbsnp_tbi)
    path(dbsnp)
    path(bed)
    path(bwa)

    output:
    tuple val(meta), path("*.svaba.unfiltered.vcf"), emit: vcf
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.patient}"
    def bwa_index = bwa ? "cp -s ${bwa}/* ."         : ""
    """
    rm ${fasta} ${fai}
    
    ${bwa_index}

    svaba run \\
        -t ${tumour_bam} \\
        -n ${normal_bam} \\
        --reference-genome ${fasta} \\
        -a ${prefix} \\
        -D ${dbsnp} \\
        -k ${bed} \\
        -p $task.cpus \\
        $args

    awk 'BEGIN {FS=OFS=\"\\\\t\"}  /^#/ {print}' ${prefix}.svaba.unfiltered.somatic.sv.vcf > ${prefix}.svaba.unfiltered.vcf

    awk 'BEGIN {FS=OFS=\"\\\\t\"}  \$1 ~ /^(1?[0-9]|2[0-2]|X|Y)\$/ {print}' ${prefix}.svaba.unfiltered.somatic.sv.vcf  >> ${prefix}.svaba.unfiltered.vcf
    
    awk 'BEGIN {FS=OFS=\"\\\\t\"}  \$1 ~ /^(1?[0-9]|2[0-2]|X|Y)\$/ {print}' ${prefix}.svaba.unfiltered.germline.sv.vcf >> ${prefix}.svaba.unfiltered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svaba: \$(svaba --version |& sed '1!d ; s/svaba //')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.svaba.unfiltered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svaba: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
