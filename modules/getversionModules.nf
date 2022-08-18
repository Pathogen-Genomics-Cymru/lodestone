// modules for the getversion workflow

process getversion {

    label 'getversion'

    publishDir "${params.output_dir}", mode: 'copy', pattern: '*.json', overwrite: 'true'

    output:
    path("version.json", emit: getversion_json)

    script:

    """
    python3 ${baseDir}/bin/software-json.py ${params.sing_dir}
    """

    stub:

    """
    touch version.json
    """
}

