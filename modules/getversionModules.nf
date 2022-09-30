// modules for the getversion workflow

process getversion {

    label 'getversion'

    publishDir "${params.output_dir}", mode: 'copy', pattern: '*.json', overwrite: 'true'

    input:
    path sing_dir

    output:
    path("version.json", emit: getversion_json)

    script:

    """
    python3 /nextflow-bin/software-json.py ${sing_dir}
    """

    stub:

    """
    touch version.json
    """
}

