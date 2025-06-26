process SERACARE_CHECKUP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'biocontainers/ubuntu:20.04' }"

    input:
    tuple val(meta), path(sv_file)

    output:
    tuple val(meta), path("*.checkup_report.txt"), emit: report
    tuple val(meta), path("*.checkup_summary.json"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    #!/bin/bash
    
    # Define the expected fusion events and their coordinates
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

    # Initialize output files
    report_file="${prefix}.checkup_report.txt"
    summary_file="${prefix}.checkup_summary.json"
    
    echo "SeraCare Quality Control Report" > \$report_file
    echo "Generated: \$(date)" >> \$report_file
    echo "Input file: ${sv_file}" >> \$report_file
    echo "=======================================" >> \$report_file
    echo "" >> \$report_file

    # Initialize JSON summary
    echo "{" > \$summary_file
    echo "  \\"sample_id\\": \\"${prefix}\\"," >> \$summary_file
    echo "  \\"input_file\\": \\"${sv_file}\\"," >> \$summary_file
    echo "  \\"timestamp\\": \\"\$(date -Iseconds)\\"," >> \$summary_file
    echo "  \\"events\\": {" >> \$summary_file

    all_passed=true
    event_count=0
    total_events=\${#value_pairs[@]}

    for key in "\${!value_pairs[@]}"; do
        value1=\$key
        value2=\${value_pairs[\$key]}
        event_name=\${event_names[\$key]}
        
        # Count occurrences of the fusion event (either direction)
        count=\$(grep -E "\$value1.*\$value2|\$value2.*\$value1" "${sv_file}" | wc -l)
        
        # Determine result status
        if [ \$count -ge 2 ]; then
            result="PASS"
            status_code=2
            echo "\$event_name event found \$count times - PASS (2)" >> \$report_file
        elif [ \$count -ge 1 ]; then
            result="WARNING"
            status_code=1
            echo "\$event_name event found \$count times - WARNING (1)" >> \$report_file
        else
            result="FAIL"
            status_code=0
            echo "\$event_name event found \$count times - FAIL (0)" >> \$report_file
            all_passed=false
        fi
        
        # Add to JSON summary
        event_count=\$((event_count + 1))
        echo "    \\"${event_name}\\": {" >> \$summary_file
        echo "      \\"coordinates\\": [\\"${value1}\\", \\"${value2}\\"]," >> \$summary_file
        echo "      \\"count\\": \$count," >> \$summary_file
        echo "      \\"status\\": \\"\$result\\"," >> \$summary_file
        echo "      \\"status_code\\": \$status_code" >> \$summary_file
        
        if [ \$event_count -lt \$total_events ]; then
            echo "    }," >> \$summary_file
        else
            echo "    }" >> \$summary_file
        fi
    done

    echo "  }," >> \$summary_file

    # Final assessment
    echo "" >> \$report_file
    echo "=======================================" >> \$report_file
    if \$all_passed; then
        overall_status="PASS"
        echo "OVERALL RESULT: PASS" >> \$report_file
        echo "All expected fusion events detected. SeraCare control sample passed QC." >> \$report_file
        echo "Pipeline ready for finalization and cleanup." >> \$report_file
    else
        overall_status="FAIL"
        echo "OVERALL RESULT: FAIL" >> \$report_file
        echo "Some expected fusion events were not detected adequately." >> \$report_file
        echo "Consider re-running the pipeline or investigating the sample quality." >> \$report_file
    fi

    # Complete JSON summary
    echo "  \\"overall_status\\": \\"\$overall_status\\"," >> \$summary_file
    echo "  \\"all_events_passed\\": \$all_passed" >> \$summary_file
    echo "}" >> \$summary_file

    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | cut -d' ' -f4)
        grep: \$(grep --version | head -n1 | cut -d' ' -f4)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.checkup_report.txt
    touch ${prefix}.checkup_summary.json
    touch versions.yml
    """
}