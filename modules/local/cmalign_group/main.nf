process CMALIGN_GROUP {
    tag "$rfam_id"
    label 'process_medium'

    conda "bioconda::infernal=1.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/infernal:1.1.5--pl5321h031d066_3' :
        'quay.io/biocontainers/infernal:1.1.5--pl5321h031d066_3' }"

    input:
    tuple val(rfam_id), path(fasta_files)
    path cm_index

    output:
    tuple val(rfam_id), path("${rfam_id}.msa.sto"), emit: stockholm
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Combine all FASTA files into one
    cat ${fasta_files} > combined.fasta

    # Extract the specific CM for this Rfam family
    cmfetch Rfam.cm ${rfam_id} > ${rfam_id}.cm

    # Align all sequences together to get a proper MSA
    cmalign \\
        --outformat Stockholm \\
        -o ${rfam_id}.msa.sto \\
        --cpu ${task.cpus} \\
        ${args} \\
        ${rfam_id}.cm \\
        combined.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        infernal: \$(cmalign -h | grep -E '^# INFERNAL' | sed 's/# INFERNAL //' | sed 's/ .*//')
    END_VERSIONS
    """
}
