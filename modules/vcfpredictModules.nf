// modules for the vcfpredict workflow

process vcfmix {

    tag {sample_name}
    label 'vcfmix'
    label 'normal_cpu'
    label 'low_memory'

    publishDir "${params.output_dir}/${sample_name}/output_vcfs", mode: 'copy', pattern: '*.json', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}/output_vcfs", mode: 'copy', pattern: '*.csv', overwrite: 'true'

    input:
    tuple val(sample_name), path(vcf)

    output:
    tuple val(sample_name), path("${sample_name}_f-stats.json"), path("${sample_name}_vcfmix-regions.csv"), emit: vcfmix_json_csv

    script:
    bcftools_vcf = "${sample_name}.bcftools.vcf"

    """
    python3 ${baseDir}/bin/vcfmix.py ${bcftools_vcf}
    """

    stub:
    vcfmix_json = "${sample_name}_f-stats.json"
    vcfmix_csv = "${sample_name}_vcfmix-regions.csv"

    """
    touch ${vcfmix_json}
    touch ${vcfmix_csv}
    """
}

process gnomon {

    tag {sample_name}
    label 'vcfmix'
    label 'normal_cpu'
    label 'low_memory'

    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*.json', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*.csv', overwrite: 'true'
    publishDir "${params.output_dir}/${sample_name}/antibiogram", mode: 'copy', pattern: '*.fasta', overwrite: 'true'

    input:
    tuple val(sample_name), path(vcf), val(isSampleTB)

    when:
    isSampleTB =~ /CREATE\_ANTIBIOGRAM\_${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}.gnomon.json"), path("${sample_name}-effects.csv"), path("${sample_name}-mutations.csv"), emit: gnomon_json_csv
    tuple val(sample_name), path("${sample_name}-fixed.fasta"), emit: gnomon_fasta

    script:
    minos_vcf = "${sample_name}.minos.vcf"

    """
    gnomon --genome_object /resources/H37rV_v3.gbk --catalogue /tuberculosis_amr_catalogues/catalogues/NC_000962.3/NC_000962.3_CRyPTIC_v1.311_GARC1_RUS.csv --vcf-file ${minos_vcf} --output_dir . --json --fasta fixed

    mv gnomon-out.json ${sample_name}.gnomon.json
    mv effects.csv ${sample_name}-effects.csv
    mv mutations.csv ${sample_name}-mutations.csv
    """

    stub:
    gnomon_json = "${sample_name}.gnomon.json"
    gnomon_fasta = "${sample_name}-fixed.fasta"
    gnomon_effects = "${sample_name}-effects.csv"
    gnomon_mutations = "${sample_name}-mutations.csv"

    """
    touch ${gnomon_json}
    touch ${gnomon_fasta}
    touch ${gnomon_effects}
    touch ${gnomon_mutations}
    """
}
