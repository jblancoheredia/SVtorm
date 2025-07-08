process DELLY {
    tag "$meta.patient"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/delly:1.2.6' :
        'blancojmskcc/delly:1.2.6' }"

    input:
    tuple val(meta) , path(tumour_bam), path(tumour_bai), path(normal_bam), path(normal_bai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fai)
    path(exclude_bed)

    output:
    tuple val(meta), path("*.{csi,tbi}")            , emit: csi
    tuple val(meta), path("*.{bcf,delly.vcf.gz}")   , emit: bcf
    tuple val(meta), path("*.delly.unfiltered.vcf") , emit: vcf
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    def suffix = task.ext.suffix ?: "bcf"
    def exclude = exclude_bed ? "--exclude ${exclude_bed}" : ""
    def bcf_output = suffix == "bcf" ? "--outfile ${prefix}.bcf" : ""
    def bcf_filter = suffix == "bcf" ? "--outfile ${prefix}.filtered.bcf" : ""
    """
    normal_name=\$(samtools view -H ${normal_bam} | grep '^@RG' | sed -n 's/.*SM:\\([^[:space:]]*\\).*/\\1/p' | head -n1)
    tumour_name=\$(samtools view -H ${tumour_bam} | grep '^@RG' | sed -n 's/.*SM:\\([^[:space:]]*\\).*/\\1/p' | head -n1)

    echo -e "\${tumour_name}\\ttumor\\n\${normal_name}\\tcontrol" > sample_file.tsv

    delly \\
        call \\
        ${args} \\
        ${bcf_output} \\
        --genome ${fasta} \\
        ${exclude} \\
        ${tumour_bam} \\
        ${normal_bam}

    delly filter \\
        -t \\
        -v 2 \\
        -m 30  \\
        -a 0.005 \\
        -f somatic \\
        -s sample_file.tsv \\
        ${bcf_filter} \\
        ${prefix}.bcf

    bcftools convert -O v -o ${prefix}.delly.vcf ${prefix}.filtered.bcf

    grep -v "^#" ${prefix}.delly.vcf | awk 'BEGIN {FS=OFS="\\t"}  \$1 ~ /^(1?[0-9]|2[0-2]|X|Y)\$/ {print}' > tmp.delly.unfiltered.vcf

    grep "^#" ${prefix}.delly.vcf > ${prefix}.delly.unfiltered.vcf

    awk 'BEGIN {FS=OFS="\\t"} /^#/ {print} (!/^#/ && \$7=="LowQual" && \$6 >= 250) 
        { \$7="PASS"; print; next } (!/^#/) { print }' tmp.delly.unfiltered.vcf >> ${prefix}.delly.unfiltered.vcf
        
    bgzip ${prefix}.delly.vcf

    rm tmp.delly.unfiltered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${patient}"
    """
    touch ${prefix}.delly.vcf.gz
    touch ${prefix}.delly.unfiltered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
    END_VERSIONS
    """
}
