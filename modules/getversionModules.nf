// modules for the getversion workflow

process getversion {

    label 'run_local'

    publishDir "${params.output_dir}", mode: 'copy', pattern: '*.json', overwrite: 'true'

    output:
    path("version.json", emit: getversion_json)

    script:

    """
    software-json.py ${params.sing_dir} ${params.config_dir}
    """

    stub:

    """
    touch version.json
    """
}

