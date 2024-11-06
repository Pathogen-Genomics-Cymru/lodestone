// modules for the getversion workflow

process getversion {

    label 'getversion'

    publishDir "${params.output_dir}", mode: 'copy', pattern: '*.json', overwrite: 'true'

    output:
    path("version.json", emit: getversion_json)

    script:

    """
    software-json.py ${params.sing_dir} ${params.config_file}
    """

    stub:

    """
    touch version.json
    """
}

