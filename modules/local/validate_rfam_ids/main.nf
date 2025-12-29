process VALIDATE_RFAM_IDS {
    tag "validate_${rfam_ids.size()}_families"
    label 'process_single'

    conda "bioconda::infernal=1.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/infernal:1.1.5--pl5321h031d066_3' :
        'quay.io/biocontainers/infernal:1.1.5--pl5321h031d066_3' }"

    input:
    val rfam_ids
    path cm_index

    output:
    path "valid_rfam_ids.txt"  , emit: valid_ids
    path "invalid_rfam_ids.txt", emit: invalid_ids
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ids_list = rfam_ids.join('\n')
    """
    # Write all Rfam IDs to check
    cat <<'EOF' > rfam_ids_to_check.txt
${ids_list}
EOF

    # Check each Rfam ID against the CM database
    touch valid_rfam_ids.txt
    touch invalid_rfam_ids.txt

    while read -r rfam_id; do
        if [ -n "\$rfam_id" ]; then
            # Try to fetch the CM - if it succeeds, the ID is valid
            if cmfetch Rfam.cm "\$rfam_id" > /dev/null 2>&1; then
                echo "\$rfam_id" >> valid_rfam_ids.txt
            else
                echo "\$rfam_id" >> invalid_rfam_ids.txt
            fi
        fi
    done < rfam_ids_to_check.txt

    # Report counts
    valid_count=\$(wc -l < valid_rfam_ids.txt)
    invalid_count=\$(wc -l < invalid_rfam_ids.txt)
    echo "Validated Rfam IDs: \$valid_count valid, \$invalid_count invalid"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        infernal: \$(cmfetch -h | grep -E '^# INFERNAL' | sed 's/# INFERNAL //' | sed 's/ .*//')
    END_VERSIONS
    """
}
