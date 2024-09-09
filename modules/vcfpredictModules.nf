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
    tag { sample_name }
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
    tag { sample_name }
    label 'medium_memory'
    label 'medium_cpu'
    label 'tbprofiler'
    
    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*.tbprofiler-out.json', overwrite: 'true'
    publishDir "${params.output_dir}${sample_name}", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    val(sample_name)
    path(minos_vcf)
    path(report_json)
    val(isSampleTB)

    output:
    tuple val(sample_name), path("${sample_name}.tbprofiler-out.json"), path("${sample_name}_report.json"), emit: tbprofiler_json
    path("${sample_name}/${sample_name}.results.json"), emit: collate_json

    when:
    isSampleTB =~ /CREATE\_ANTIBIOGRAM\_${sample_name}/

    script:
    error_log = "${sample_name}_err.json"
    tbprofiler_json = "${sample_name}.tbprofiler-out.json"
    
    """
    bgzip ${minos_vcf}
    
    mkdir tmp
    tb-profiler profile --vcf ${minos_vcf}.gz --threads ${task.cpus} --temp tmp --prefix ${sample_name}
    
    mv results ${sample_name}
    cp ${sample_name}/${sample_name}.results.json ${tbprofiler_json}
    
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

process ntmprofiler {
    tag {sample_name}
    label 'low_memory'
    label 'low_cpu'
    label 'ntmprofiler'
   
    input:
    tuple val(sample_name), path(fq1), path(fq2), path(report_json), val(isSampleTB)
    
    output:
    path("${sample_name}.ntmprofiler-out.json"), path("${sample_name}_report.json"), emit: ntmprofiler_json
    path("${sample_name}.results.json"), emit: collate_json

    when:
    isSampleTB != /CREATE\_ANTIBIOGRAM\_${sample_name}/

    script:
    error_log = "${sample_name}_err.json"
    ntmprofiler_json = "${sample_name}.ntmprofiler-out.json"

    """
    mkdir tmp
    ntm-profiler profile -1 $fq1 -2 $fq2 --threads ${task.cpus} --temp tmp --prefix ${sample_name}
    
    cp ${sample_name}.results.josn ${ntmprofiler_json}

    cp ${sample_name}_report.json ${sample_name}_report_previous.json

    echo '{"complete":"workflow complete without error"}' | jq '.' > ${error_log}

    jq -s ".[0] * .[1] * .[2]" ${error_log} ${sample_name}_report_previous.json  ${ntmprofiler_json} > ${report_json}
    """
}

process tbtamr {
    tag { sample_name }
    label 'medium_memory'
    label 'medium_cpu'
    label 'tbtamr'
    
    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*.tbtamr-out.json', overwrite: 'true'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), path(report_json), val(isSampleTB)

    output:
    tuple val(sample_name), path("${sample_name}.tbtamr-out.json"), path("${sample_name}_report.json"), emit: tbtamr_json
    path(sample_name), emit: collate_json

    when:
    isSampleTB =~ /CREATE\_ANTIBIOGRAM\_${sample_name}/

    script:
    error_log = "${sample_name}_err.json"
    tbtamr_json = "${sample_name}.tbtamr-out.json"
    
    """
    tbtamr run -r1 $fq1 -r2 $fq2 --prefix ${sample_name}
    
    cp ${sample_name}/tbtamr.json ${tbtamr_json}
    
    cp ${sample_name}_report.json ${sample_name}_report_previous.json

    echo '{"complete":"workflow complete without error"}' | jq '.' > ${error_log}

    #tidy up report so we can combine
    cp ${tbtamr_json} ${sample_name}_tbtamr.json
    sed -i '1d;\$d' ${tbtamr_json}
    sed -i 's/Seq_ID/resistance_profiler/g' ${tbtamr_json}

    jq -s ".[0] * .[1] * .[2]" ${error_log} ${sample_name}_report_previous.json  ${tbtamr_json} > ${report_json}
    """

    stub:
    """
    touch ${sample_name}.tbtamr-out.json
    touch ${sample_name}_report.json
    """
}

process tbtamr_collate{
    label 'medium_memory'
    label 'medium_cpu'
    label 'tbtamr'

    publishDir "${params.output_dir}", mode: 'copy', overwrite: 'true', pattern: 'tbtamr.csv', stageAs: "tbtamr.variants.csv"

    input:
    path(files), stageAs: "results/*"

    output:
    path("tbtamr.csv")

    script:
    """
    tbtamr collate
    """
}

process tbprofiler_collate{
    label 'medium_memory'
    label 'medium_cpu'
    label 'tbprofiler'

    publishDir "${params.output_dir}", mode: 'copy', overwrite: 'true', pattern: 'tbprofiler.variants.csv'

    input:
    path(files), stageAs: "results/*"
    
    output:
    path("tbprofiler.variants.csv")

    script:
    """
    tb-profiler collate
    """
}

process ntmprofiler_collate{
    label 'medium_memory'
    label 'medium_cpu'
    label 'ntmprofiler'

    publishDir "${params.output_dir}", mode: 'copy', owerwrite: 'true', pattern: ' ntmprofiler.collate.txt.variants.csv'

    input:
    path(files)
    
    output:
    path('ntmprofiler.collate.txt.variants.csv')

    script:
    """
    ntm-profiler collate 
    """
}

process add_allelic_depth {
    tag { sample_name }
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
