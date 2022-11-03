// modules for the getversion workflow
params.output_dir = "${params.output_dir}/${workflow.sessionId}"

process getversion {

    label 'getversion'

    publishDir "${params.output_dir}/${workflow.sessionId}", mode: 'copy', pattern: '*.json', overwrite: 'true'

    input:
    path sing_dir

    output:
    path("version.json", emit: getversion_json)

    script:

    """
    software-json.py ${sing_dir}
    """

    stub:

    """
    touch version.json
    """
}

