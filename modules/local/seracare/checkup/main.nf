process SERACARE_CHECKUP {
    tag "$meta.patient"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'biocontainers/ubuntu:20.04' }"

    input:
    tuple val(meta), path(sv_file)

    output:
    tuple val(meta), path("*_SeraCareCheckUpReport.txt"), emit: report
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    #!/bin/bash
    
    declare -A value_pairs=(
        ["294472"]="425227"
        ["149784"]="117647"
        ["515842"]="436105"
    )
    
    declare -A event_names=(
        ["294472"]="ALK-EML4"
        ["149784"]="CD74-ROS1"
        ["515842"]="RET-NCOA4"
    )
    
    report_file="${prefix}_SeraCareCheckUpReport.txt"
    
    echo "SeraCare Quality Control Report" > \$report_file
    echo "Generated: \$(date)" >> \$report_file
    echo "Input file: ${sv_file}" >> \$report_file
    echo "====================================================================================================" >> \$report_file
    echo "" >> \$report_file
    
    all_passed=true
    for key in "\${!value_pairs[@]}"; do
        value1=\$key
        value2=\${value_pairs[\$key]}
        event_name=\${event_names[\$key]}
        count=\$(grep -E "\$value1.*\$value2|\$value2.*\$value1" "\$file" | wc -l)
        if [ \$count -ge 2 ]; then
            result=2
            echo "\$event_name event found \$count times - PASS (2)" >> \$report_file
        elif [ \$count -ge 1 ]; then
            result=1
            echo "\$event_name event found \$count times - WARNING (1)" >> \$report_file
        else
            result=0
            echo "\$event_name event found \$count times - FAIL (0)" >> \$report_file
            all_passed=false
        fi
    done
    if \$all_passed; then
        echo "All pairs passed in file \$file. Ready for finalizing and cleanup." >> \$report_file
    else
        echo "Some pairs did not pass in file \$file. Consider re-running the pipeline." >> \$report_file
    fi
    echo "----------------------------------------------------------------------------------------------------"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | cut -d' ' -f4)
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient}"
    """
    touch ${prefix}_SeraCareCheckUpReport.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | cut -d' ' -f4)
    END_VERSIONS
    """
}