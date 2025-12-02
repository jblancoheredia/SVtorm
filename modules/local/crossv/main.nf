process CROSSV {
    tag "$meta.patient"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/crossv_svtorm:1.0.3':
        'blancojmskcc/crossv_svtorm:1.0.3' }"

    input:
    tuple val(meta), path(tumour_bam)  , 
                     path(tumour_bai)  ,
                     path(normal_bam)  , 
                     path(normal_bai)  ,
                     path(delly_vcf)   ,
                     path(gridss_vcf)  ,
                     path(manta_vcf)   ,
                     path(recall_vcf)  ,
                     path(svaba_vcf)   ,
                     path(tsv)         ,
                     path(bed_file)

    output:
    tuple val(meta), path("*_ANNOTE_SV_INN.tsv"), emit: annote
    tuple val(meta), path("*_CrosSV.json")      , emit: json
    tuple val(meta), path("*_CrosSV.bam")       , emit: bam,  optional: true
    tuple val(meta), path("*_CrosSV.bai")       , emit: bai,  optional: true
    tuple val(meta), path("*_CrosSV.log")       , emit: log
    tuple val(meta), path("*_CrosSV.tsv")       , emit: tsv
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bed = task.ext.bed ?: bed_file
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    ln -sf \$(basename ${delly_vcf}) ${prefix}_DELLY__SV_UNF.vcf
    ln -sf \$(basename ${gridss_vcf}) ${prefix}_GRIDSS_SV_UNF.vcf
    ln -sf \$(basename ${manta_vcf}) ${prefix}_MANTA__SV_UNF.vcf
    ln -sf \$(basename ${recall_vcf}) ${prefix}_RECALL_SV_UNF.vcf
    ln -sf \$(basename ${svaba_vcf}) ${prefix}_SVABA__SV_UNF.vcf
    
    CrosSV \\
      --slop 100 \\
      --bed ${bed} \\
      --relax-type  \\
      --final-tsv ${tsv} \\
      --bam ${tumour_bam} \\
      --bam-out-prefix ${prefix} \\
      --log ${prefix}_CrosSV.log  \\
      --out-tsv ${prefix}_CrosSV.tsv \\
      --out-json ${prefix}_CrosSV.json \\
      --vcf_list <(ls ${prefix}_*_SV_UNF.vcf)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crossv: "1.0.3"
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    touch ${prefix}_ANNOTE_SV_INN.tsv
    touch ${prefix}_CrosSV.json
    touch ${prefix}_CrosSV.bam
    touch ${prefix}_CrosSV.bai
    touch ${prefix}_CrosSV.log
    touch ${prefix}_CrosSV.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crossv: "1.0.3"
    END_VERSIONS
    """
}
