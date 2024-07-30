// modules for the vcfpredict workflow

process vcfmix {

    tag {sample_name}
    label 'vcfpredict'
    label 'low_memory'
    label 'low_cpu'

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
    run-vcfmix.py ${bcftools_vcf}

    cp ${sample_name}_report.json ${sample_name}_report_previous.json

    jq -s ".[0] * .[1]" ${sample_name}_report_previous.json ${sample_name}_f-stats.json > ${report_json}

    if [ ${params.resistance_profiler} == "none" ]; then echo '{"complete":"workflow complete without error"}' | jq '.' > ${error_log} && jq -s ".[0] * .[1] * .[2]" ${error_log} ${sample_name}_report_previous.json ${sample_name}_f-stats.json > ${report_json}; fi
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

process tbprofiler_update_db {
    label 'low_memory'
    label 'low_cpu'
    label 'tbprofiler'

    input:
    path(reference)

    script:
    """
    mkdir tmp
    tb-profiler update_tbdb --match_ref $reference --temp tmp
    """
}

process tbprofiler {
    tag {sample_name}
    label 'medium_memory'
    label 'medium_cpu'
    label 'tbprofiler'
    
    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*.tbprofiler-out.json', overwrite: 'true'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    val(sample_name)
    path(minos_vcf)
    path(report_json)
    val(isSampleTB)

    output:
    tuple val(sample_name), path("${sample_name}.tbprofiler-out.json"), path("${sample_name}_report.json"), emit: tbprofiler_json

    when:
    isSampleTB =~ /CREATE\_ANTIBIOGRAM\_${sample_name}/

    script:
    error_log = "${sample_name}_err.json"
    tbprofiler_json = "${sample_name}.tbprofiler-out.json"
    
    """
    bgzip ${minos_vcf}
    
    mkdir tmp
    tb-profiler profile --vcf ${minos_vcf}.gz --threads ${task.cpus} --temp tmp
    
    mv results/tbprofiler.results.json ${tbprofiler_json}
    
    cp ${sample_name}_report.json ${sample_name}_report_previous.json

    echo '{"complete":"workflow complete without error"}' | jq '.' > ${error_log}

    jq -s ".[0] * .[1] * .[2]" ${error_log} ${sample_name}_report_previous.json  ${tbprofiler_json} > ${report_json}
    """

    stub:
    """
    touch ${sample_name}.tbprofiler-out.json
    touch ${sample_name}_report.json
    """
}

process add_allelic_depth {
    tag {sample_name}
    label 'low_memory'
    label 'low_cpu'
    label 'tbprofiler'
    
    input:
    val(sample_name)
    path(minos_vcf)
    path(bam)
    path(reference)
    val(isSampleTB)
    
    output:
    path("${sample_name}_allelic_depth.minos.vcf")

    when:
    isSampleTB =~ /CREATE\_ANTIBIOGRAM\_${sample_name}/
    
    script:
    """
    samtools faidx $reference
    samtools dict $reference -o ${reference.baseName}.dict
    
    mkdir tmp
    
    gatk VariantAnnotator -R $reference -I $bam -V $minos_vcf -A DepthPerAlleleBySample -O ${sample_name}_allelic_depth.minos.vcf --tmp-dir tmp
    """

    stub:
    """
    touch ${sample_name}_allelic_depth.minos.vcf 
    """
    
}

process gnomonicus {

    tag {sample_name}
    label 'vcfpredict'
    label 'low_memory'
    label 'low_cpu'

    errorStrategy 'ignore'

    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*.gnomonicus-out.json', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*.csv', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*.fasta', overwrite: 'true'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(vcf), val(isSampleTB), path(report_json)
    path(genbank)
    when:
    isSampleTB =~ /CREATE\_ANTIBIOGRAM\_${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}.gnomonicus-out.json"), path("${sample_name}_report.json"), emit: gnomon_json
    tuple val(sample_name), path("${sample_name}.effects.csv"), path("${sample_name}.mutations.csv"), emit: gnomon_csv optional true
    tuple val(sample_name), path("*-fixed.fasta"), emit: gnomon_fasta
    path("${sample_name}_err.json", emit: gnomon_log)
    path ("${sample_name}_report.json", emit: gnomon_report)

    script:
    minos_vcf = "${sample_name}.minos.vcf"
    error_log = "${sample_name}_err.json"

    """
    gnomonicus --genome_object ${genbank} --catalogue ${params.amr_cat} --vcf_file ${minos_vcf} --output_dir . --json --fasta fixed

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
    label 'low_memory'
    label 'low_cpu'
    
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
