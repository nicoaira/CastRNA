process GENERATE_MSA {
    tag "$rfam_id"
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(rfam_id), path(stockholm)

    output:
    tuple val(rfam_id), path("${rfam_id}.msa.sto")  , emit: stockholm
    tuple val(rfam_id), path("${rfam_id}.aligned.fa"), emit: fasta
    tuple val(rfam_id), path("${rfam_id}.structures.tsv"), emit: tsv
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate_msa_outputs.py \\
        --input ${stockholm} \\
        --output-sto ${rfam_id}.msa.sto \\
        --output-fasta ${rfam_id}.aligned.fa \\
        --output-tsv ${rfam_id}.structures.tsv \\
        --rfam-id ${rfam_id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
