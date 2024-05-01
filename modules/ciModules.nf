process formatInput {
    label 'run_local'

    input:
    tuple val(sample_name), path(fq1), path(fq2)

    output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, path("software.json"), emit: inputfqs

    script:
    """
    touch software.json
    echo /${sample_name}/
    """
}