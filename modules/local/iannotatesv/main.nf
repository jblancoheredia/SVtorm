process IANNOTATESV {
    tag "$meta.patient"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/iannotatesv:1.2.1':
        'blancojmskcc/iannotatesv:1.2.1' }"

    input:
    tuple val(meta),  path(filtered_vcf), path(filtered_vcf_index), path(filtered_tsv), path(annote_input)

    output:
    tuple val(meta), file("*_Annotated.txt")    , emit: txt
    tuple val(meta), file("*_Annotated.json")   , emit: json
    tuple val(meta), file("*_Annotated.xlsx")   , emit: xlsx
    tuple val(meta), file("*_functional.txt")   , emit: ftxt
    tuple val(meta), file("*_SOMTIC_SV_OUT.tsv"), emit: tsv
    tuple val(meta), file("*_SOMTIC_SV_ANN.tsv"), emit: ann
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    python /opt/iAnnotateSV/iAnnotateSV/iAnnotateSV.py \\
    -i ${annote_input} \\
    -ofp ${prefix} \
    -o . \
    -r hg19 \
    -d 3000

    paste ${filtered_tsv} ${prefix}_Annotated.txt > ${prefix}_SOMTIC_SV_ANN.tsv

    innput_file=\"${prefix}_SOMTIC_SV_ANN.tsv\"
    output_file=\"${prefix}_SOMTIC_SV_OUT.tsv\"

    supp_col=17
    gene1col=26
    gene2col=29

    int_threshold=3

    declare -a allow_list=(\"ALK\" \"BRAF\" \"EGFR\" \"ETV6\" \"FGFR2\" \"FGFR3\" \"MET\" \"NTRK1\" \"RET\" \"ROS1\")

    allow_pattern=\$(IFS=\"|\"; echo \"\${allow_list[*]}\")

    awk -v x_col=\"\$supp_col\" -v y_col=\"\$gene1col\" -v z_col=\"\$gene2col\" -v threshold=\"\$int_threshold\" -v pattern=\"\$allow_pattern\" '
    {FS=\"\\\\t\"}{
      x_value = \$x_col;
      y_value = \$y_col;
      z_value = \$z_col;

      print \"Processing line: \" NR;
      print \"x_value: \" x_value, \"y_value: \" y_value, \"z_value: \" z_value;
      print \"Matches pattern? y_value: \" (y_value ~ pattern), \"z_value: \" (z_value ~ pattern);

      if (x_value > threshold || y_value ~ pattern || z_value ~ pattern) {
        print \"Line \" NR \" passed the filter.\";
        print \$0 > \"'\$output_file'\"
    } else {
        print \"Line \" NR \" did not pass the filter.\";
    }
    }' \"\$innput_file\"

    echo \"Filtering complete. Results saved to \$output_file.\"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iannotatesv: "1.2.1"
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    touch ${prefix}_SOMTIC_SV_OUT.tsv
    touch ${prefix}_SOMTIC_SV_ANN.tsv
    touch ${prefix}_Annotated.txt
    touch ${prefix}_Annotated.json
    touch ${prefix}_Annotated.xlsx
    touch ${prefix}_functional.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iannotatesv: "1.2.1"
    END_VERSIONS
    """
}
