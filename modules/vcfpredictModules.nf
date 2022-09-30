// modules for the vcfpredict workflow

process vcfmix {

    tag {sample_name}
    label 'vcfpredict'
    label 'normal_cpu'
    label 'low_memory'

    publishDir "${params.output_dir}/${sample_name}/output_vcfs", mode: 'copy', pattern: '*.json', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}/output_vcfs", mode: 'copy', pattern: '*.csv', overwrite: 'true'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*.err'

    input:
    tuple val(sample_name), path(vcf)

    output:
    tuple val(sample_name), path("${sample_name}_f-stats.json"), path("${sample_name}_vcfmix-regions.csv"), emit: vcfmix_json_csv
    path "${sample_name}.err", emit: vcfmix_log optional true

    script:
    bcftools_vcf = "${sample_name}.bcftools.vcf"
    error_log = "${sample_name}.err"

    """
    python3 /nextflow-bin/vcfmix.py ${bcftools_vcf}

    if [ ${params.gnomon} == "no" ]; then printf "workflow complete without error" >> ${error_log}; fi
    """

    stub:
    vcfmix_json = "${sample_name}_f-stats.json"
    vcfmix_csv = "${sample_name}_vcfmix-regions.csv"
    error_log = "${sample_name}.err"

    """
    touch ${vcfmix_json}
    touch ${vcfmix_csv}
    touch ${error_log}
    """
}

process gnomon {

    tag {sample_name}
    label 'vcfpredict'
    label 'normal_cpu'
    label 'low_memory'

    errorStrategy 'ignore'

    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*.json', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*.csv', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*.fasta', overwrite: 'true'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*.err'

    input:
    tuple val(sample_name), path(vcf), val(isSampleTB)
    path resources_dir

    when:
    isSampleTB =~ /CREATE\_ANTIBIOGRAM\_${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}.gnomon-out.json"), path("${sample_name}.effects.csv"), path("${sample_name}.mutations.csv"), emit: gnomon_json_csv
    tuple val(sample_name), path("*-fixed.fasta"), emit: gnomon_fasta
    path("${sample_name}.err", emit: gnomon_log)

    script:
    minos_vcf = "${sample_name}.minos.vcf"
    error_log = "${sample_name}.err"

    """
    gnomon --genome_object ${resources_dir}/H37rV_v3.gbk --catalogue ${params.amr_cat} --vcf_file ${minos_vcf} --output_dir . --json --fasta fixed

    printf "workflow complete without error" >> ${error_log} 
    """

    stub:
    gnomon_json = "${sample_name}.gnomon-out.json"
    gnomon_fasta = "${sample_name}-fixed.fasta"
    gnomon_effects = "${sample_name}.effects.csv"
    gnomon_mutations = "${sample_name}.mutations.csv"
    error_log = "${sample_name}.err"

    """
    touch ${gnomon_json}
    touch ${gnomon_fasta}
    touch ${gnomon_effects}
    touch ${gnomon_mutations}
    touch ${error_log}
    """
}
