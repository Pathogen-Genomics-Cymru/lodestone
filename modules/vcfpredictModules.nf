// modules for the vcfpredict workflow

process vcfmix {

    tag {sample_name}
    label 'vcfpredict'
    label 'normal_cpu'
    label 'low_memory'

    errorStrategy 'ignore'

    publishDir "${params.output_dir}/${sample_name}/output_vcfs", mode: 'copy', pattern: '*_f-stats.json', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}/output_vcfs", mode: 'copy', pattern: '*.csv', overwrite: 'true'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(vcf), path(report_json)

    output:
    tuple val(sample_name), path("${sample_name}_f-stats.json"), emit: vcfmix_json
    tuple val(sample_name), path("${sample_name}_f-stats.json"), path("${sample_name}_vcfmix-regions.csv"), emit: vcfmix_json_csv
    path "${sample_name}_err.json", emit: vcfmix_log optional true
    path ("${sample_name}_report.json", emit: vcfmix_report)

    script:
    bcftools_vcf = "${sample_name}.bcftools.vcf"
    error_log = "${sample_name}_err.json"

    """
    python3 ${baseDir}/bin/vcfmix.py ${bcftools_vcf}

    cp ${sample_name}_report.json ${sample_name}_report_previous.json

    jq -s ".[0] * .[1]" ${sample_name}_report_previous.json ${sample_name}_f-stats.json > ${report_json}

    if [ ${params.gnomonicus} == "no" ]; then echo '{"complete":"workflow complete without error"}' | jq '.' > ${error_log} && jq -s ".[0] * .[1] * .[2]" ${error_log} ${sample_name}_report_previous.json ${sample_name}_f-stats.json > ${report_json}; fi
    """

    stub:
    vcfmix_json = "${sample_name}_f-stats.json"
    vcfmix_csv = "${sample_name}_vcfmix-regions.csv"
    error_log = "${sample_name}_err.json"

    """
    touch ${vcfmix_json}
    touch ${vcfmix_csv}
    touch ${error_log}
    """
}

process gnomonicus {

    tag {sample_name}
    label 'vcfpredict'
    label 'normal_cpu'
    label 'low_memory'

    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*.gnomonicus-out.json', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*.csv', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*.fasta', overwrite: 'true'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(vcf), val(isSampleTB), path(report_json)

    when:
    isSampleTB =~ /CREATE\_ANTIBIOGRAM\_${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}.gnomonicus-out.json"), path("${sample_name}_report.json"), emit: gnomon_json
    tuple val(sample_name), path("${sample_name}.gnomonicus-out.json"), path("${sample_name}.effects.csv"), path("${sample_name}.mutations.csv"), emit: gnomon_json_csv
    tuple val(sample_name), path("*-fixed.fasta"), emit: gnomon_fasta
    path("${sample_name}_err.json", emit: gnomon_log)
    path ("${sample_name}_report.json", emit: gnomon_report)

    script:
    minos_vcf = "${sample_name}.minos.vcf"
    error_log = "${sample_name}_err.json"

    """
    gnomonicus --genome_object ${baseDir}/resources/H37rV_v3.gbk --catalogue ${params.amr_cat} --vcf_file ${minos_vcf} --output_dir . --json --fasta fixed

    cp ${sample_name}_report.json ${sample_name}_report_previous.json

    echo '{"complete":"workflow complete without error"}' | jq '.' > ${error_log}

    jq -s ".[0] * .[1] * .[2]" ${error_log} ${sample_name}_report_previous.json ${sample_name}.gnomonicus-out.json > ${report_json}
    """

    stub:
    gnomonicus_json = "${sample_name}.gnomonicus-out.json"
    gnomonicus_fasta = "${sample_name}-fixed.fasta"
    gnomonicus_effects = "${sample_name}.effects.csv"
    gnomonicus_mutations = "${sample_name}.mutations.csv"
    error_log = "${sample_name}_err.json"

    """
    touch ${gnomonicus_json}
    touch ${gnomonicus_fasta}
    touch ${gnomonicus_effects}
    touch ${gnomonicus_mutations}
    touch ${error_log}
    """
}

process finalJson {

    tag {sample_name}
    label 'vcfpredict'
    label 'normal_cpu'
    label 'low_memory'

    errorStrategy 'ignore'

    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*_report.json'

    input:
    tuple val(sample_name), path(vcfmix_json), path(gnomon_json), path(report_json)

    output:
    tuple val(sample_name), path("${sample_name}_report.json"), emit: final_json

    script:
    """
    cp ${sample_name}_report.json ${sample_name}_report_previous.json

    jq -s ".[0] * .[1]" ${sample_name}_report_previous.json ${vcfmix_json} > ${report_json}
    """

    stub:
    report_json = "${sample_name}_report.json"

    """
    touch ${report_json}
    """

}
