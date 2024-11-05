process identifyBacterialContaminants {
    /**
    * @QCcheckpoint if urllist.txt is empty, there are no contaminant genomes to download, so skip next process
    */


    tag "${sample_name } pass ${pass}"
    label 'preprocessing'
    label 'normal_cpu'
    label 'medium_memory'

    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP", mode: 'copy', pattern: '*_species_in_samp*.json'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), path(software_json), path(afanc_json), val(enough_myco_reads), path(kraken_report), path(kraken_json)
    val(resources)
    path(refseq)
    val(pass)

    when:
    enough_myco_reads =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_urllist.txt"), stdout, path(software_json), path("${sample_name}_species_in_sample_previous.json"), emit: contam_list optional true
    tuple val(sample_name), path("${sample_name}_species_in_sample_previous.json"), stdout, path(software_json), emit: prev_sample_json optional true
    tuple val(sample_name), path("${sample_name}_species_in_sample_pass_${pass}.json"), stdout, emit: sample_json
    tuple val(sample_name), path("${sample_name}_nocontam_1.fq.gz"), path("${sample_name}_nocontam_2.fq.gz"), path(software_json), path("${sample_name}_species_in_sample_pass_${pass}.json"), stdout, emit: nocontam_fqs optional true
    path "${sample_name}_err.json", emit: contam_log optional true
    path "${sample_name}_report.json", emit: contam_report optional true

    script:
    error_log = "${sample_name}_err.json"
    report_json = "${sample_name}_report.json"
    report_pass_json = "${sample_name}_species_in_sample_pass_${pass}.json"

    """
    identify_tophit_and_contaminants2.py ${afanc_json} ${kraken_json} ${refseq} ${params.species} ${params.unmix_myco} ${resources} null ${params.permissive} ${pass}

    contam_to_remove=\$(jq -r '.summary_questions.are_there_contaminants' ${report_pass_json})
    acceptable_species=\$(jq -r '.summary_questions.is_the_top_species_appropriate' ${report_pass_json})
    top_hit=\$(jq -r '.top_hit.name' ${report_pass_json})

    if [ \$contam_to_remove == 'yes' ]; then 
        cp ${report_pass_json} ${sample_name}_species_in_sample_previous.json
    fi

    if [ \$contam_to_remove == 'yes' ]; then
        printf "NOW_DECONTAMINATE_${sample_name}"
    elif [ \$contam_to_remove == 'no' ] && [ \$acceptable_species == 'yes' ] && [ ${pass} == 1 ]; then 
        printf "NOW_ALIGN_TO_REF_${sample_name}" && mv $fq1 ${sample_name}_nocontam_1.fq.gz && mv $fq2 ${sample_name}_nocontam_2.fq.gz; 
    elif [ \$contam_to_remove == 'no' ] && [ \$acceptable_species == 'no' ]; then 
        jq -n --arg key "\$top_hit" '{"error": ("top hit " + \$key + " does not have a reference genome. Sample will not proceed beyond preprocessing workflow.")}' > ${error_log} && \
        jq -s ".[0] * .[1] * .[2]" ${software_json} ${error_log} ${report_pass_json} > ${report_json}
    fi
    """

    stub:
    error_log = "${sample_name}_err.json"

    """
    touch ${sample_name}_species_in_sample_pass_${pass}.json
    touch ${sample_name}_species_in_sample_previous.json
    touch ${sample_name}_urllist.txt
    touch ${error_log}
    printf ${params.identifyBacContam_rundecontam}
    """
}

process downloadContamGenomes {
    /**
    * @QCcheckpoint confirm that we could download every genome in the list of contaminants
    */

    tag "${sample_name } pass ${pass}"
    label 'preprocessing'
    label 'low_cpu'
    label 'medium_memory'

    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(contam_list), val(run_decontaminator), path(software_json), path(prev_species_json)
    val(pass)

    when:
    run_decontaminator =~ /NOW\_DECONTAMINATE\_${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_contaminants.fa"), stdout, emit: contam_fa
    path "${sample_name}_err.json", emit: downcontam_log optional true
    path "${sample_name}_report.json", emit: downcontam_report optional true

    script:
    contaminant_fa = "${sample_name}_contaminants.fa"
    error_log = "${sample_name}_err.json"
    report_json = "${sample_name}_report.json"

    """
    wget -i ${contam_list} --spider -nv -a linktestlog.txt 2>&1
    grep -o 'https://.*fna.gz' linktestlog.txt > confirmedurllist.txt

    mkdir contam_dir
    cd contam_dir

    wget -i ../confirmedurllist.txt

    gunzip *.gz
    cat *.fna > ../${contaminant_fa}
    rm -rf *.fna
    cd ..

    num_urls_in=\$(cat $contam_list | wc -l)
    num_urls_out=\$(cat confirmedurllist.txt | wc -l)

    rm -rf linktestlog.txt confirmedurllist.txt
    rm -r contam_dir

    if (( \$num_urls_in == \$num_urls_out )); then printf "${sample_name}"; else jq -n --arg key1 "\$num_urls_in" --arg key2 "\$num_urls_out" '{"error": ("there were " + \$key1 + " contaminant genomes but only " +  \$key2 +  " could be downloaded")}' > ${error_log} && printf "fail" && jq -s ".[0] * .[1] * .[2]" ${software_json} ${error_log} ${prev_species_json} > ${report_json}; fi
    """

    stub:
    contaminant_fa = "${sample_name}_contaminants.fa"
    error_log = "${sample_name}_err.json"

    """
    printf ${params.downloadContamGenomes_fapass}
    touch ${contaminant_fa}
    touch ${error_log}
    """
}

process mapToContamFa {
    /**
    * @QCcheckpoint none
    */

    tag "${sample_name } pass ${pass}"
    label 'preprocessing'
    label 'normal_cpu'
    label 'high_memory'

    publishDir "${params.output_dir}/$sample_name/output_reads", mode: 'copy', pattern: '*.fq.gz', overwrite: 'true'

    input:
    tuple val(sample_name), path(fq1), path(fq2), path(software_json), path(contam_fa), val(does_fa_pass)
    val(pass)    

    when:
    does_fa_pass =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_cleaned_1.fq.gz"), path("${sample_name}_cleaned_2.fq.gz"), path(software_json), emit: reClassification_fqs
    path(software_json), emit: software_json

    script:
    bam = "${sample_name}.bam"
    decontam_fq1 = "${sample_name}_cleaned_1.fq"
    decontam_fq2 = "${sample_name}_cleaned_2.fq"

    """
    echo ${task.index}
    bwa index ${contam_fa}
    bwa mem -t ${task.cpus} -M ${contam_fa} ${fq1} ${fq2} | samtools view -f 4 -f 8 -Shb - > ${bam}

    samtools fastq -1 ${decontam_fq1} -2 ${decontam_fq2} ${bam}
    rm -rf ${bam}

    gzip -f ${decontam_fq1}
    gzip -f ${decontam_fq2}
    """

    stub:
    decontam_fq1 = "${sample_name}_cleaned_1.fq"
    decontam_fq2 = "${sample_name}_cleaned_2.fq"

    """
    touch ${decontam_fq1}.gz
    touch ${decontam_fq2}.gz
    """
}

process reKraken {
    /**
    * @QCcheckpoint none
    */


    tag "${sample_name } pass ${pass}"
    label 'preprocessing'
    label 'normal_cpu'
    label 'high_memory'

    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP_and_postContamRemoval", mode: 'copy', pattern: '*_kraken*'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), path(software_json)
    path(database)
    val(pass)

    output:
    tuple val(sample_name), path("${sample_name}_kraken_report.txt"), path("${sample_name}_kraken_report.json"), emit: reKraken_report

    script:
    kraken2_report = "${sample_name}_kraken_report.txt"
    kraken2_json = "${sample_name}_kraken_report.json"
    kraken2_read_classification = "${sample_name}_read_classifications.txt"

    """
    kraken2 --threads ${task.cpus} --db . --output ${kraken2_read_classification} --report ${kraken2_report} --paired $fq1 $fq2

    parse_kraken_report2.py ${kraken2_report} ${kraken2_json}     ${params.kraken.kraken_percent_threshold} \
    ${params.kraken.kraken_n_reads_threshold} ${params.permissive}

    rm -rf ${sample_name}_read_classifications.txt
    """

    stub:
    kraken2_report = "${sample_name}_kraken_report.txt"
    kraken2_json = "${sample_name}_kraken_report.json"
    kraken2_read_classification = "${sample_name}_read_classifications.txt"

    """
    touch ${kraken2_report}
    touch ${kraken2_json}
    touch ${kraken2_read_classification}
    """
}

process reAfanc {
    /**
    * @QCcheckpoint none
    */


    


    tag "${sample_name } pass ${pass}"
    label 'preprocessing'
    label 'normal_cpu'
    label 'medium_memory'
    label 'retry_afanc'

    errorStrategy 'ignore'

    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP_and_postContamRemoval", mode: 'copy', pattern: '*_afanc*.json'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), path(software_json)
    path(afanc_myco_db)
    val(pass)

    output:
    tuple val(sample_name), path("${sample_name}_afanc_report.json"), emit: reAfanc_report
    path "${sample_name}_afanc_original.json", emit: reAfanc_original
    // tuple val(sample_name), path(fq1), path(fq2), stdout, emit: afanc_fqs

    script:
    afanc_report = "${sample_name}_afanc_report.json"

    """
    afanc screen ${afanc_myco_db} ${fq1} ${fq2}-p ${params.afanc.afanc_percent_threshold} \
        -n ${params.afanc.afanc_n_reads_threshold} -o ${sample_name} -t ${task.cpus} -v ${afanc_myco_db}/lineage_profiles/TB_variants.tsv 

    cp ${sample_name}/${sample_name}.json ${sample_name}_afanc_original.json
    reformat_afanc_json.py ${sample_name}/${sample_name}.json
    printf ${sample_name}
    """

    stub:
    afanc_report = "${sample_name}_afanc_report.json"
    afanc_original =  "${sample_name}_afanc_original.json"

    """
    touch ${afanc_report}
    touch ${afanc_original}
    printf ${sample_name}
    """
}

process reMykrobe {
    /**
    * @QCcheckpoint none
    */
    

    

    tag "${sample_name } pass ${pass}"
    label 'preprocessing'
    label 'normal_cpu'
    label 'low_memory'

    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP_and_postContamRemoval", mode: 'copy', pattern: '*_mykrobe_report.json'
    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP_and_postContamRemoval", mode: 'copy', pattern: '*_mykrobe_report.csv'

    input:
    tuple val(sample_name), path(fq1), path(fq2), path(software_json)
    val(pass)

    output:
    tuple val(sample_name), path("${sample_name}_mykrobe_report.json"), emit: reMykrobe_report

    script:
    mykrobe_report = "${sample_name}_mykrobe_report"

    """
    mykrobe predict --sample ${sample_name} --species tb --threads ${task.cpus} --format json_and_csv --output ${mykrobe_report} -1 $fq1 $fq2
    """

    stub:
    mykrobe_report = "${sample_name}_mykrobe_report.json"

    """
    touch ${mykrobe_report}
    """
}

process summarise {
    /**
    * @QCcheckpoint checks whether there are still contaminants and if top species hit is supported
    */

    

    tag "${sample_name } pass ${pass}"
    label 'preprocessing'
    label 'low_cpu'
    label 'medium_memory'

    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP_and_postContamRemoval", mode: 'copy', pattern: '*_species_in_sample.json'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(afanc_json), path(kraken_report), path(kraken_json), path(prev_species_json), val(decontam), path(software_json)
    val(resources)
    path(refseq)
    val(pass)

    output:
    tuple val(sample_name), path("${sample_name}_species_in_sample_pass_${pass}.json"), stdout, emit: summary_json
    stdout emit: do_we_break
    path "${sample_name}_err.json", emit: summary_log optional true
    path "${sample_name}_pass_${pass}_report.json", emit: summary_report optional true
    val(pass), emit: pass_number

    script:
    error_log = "${sample_name}_err.json"
    report_json = "${sample_name}_pass_${pass}_report.json"
    species_in_sample = "${sample_name}_species_in_sample_pass_${pass}.json"
    """
    identify_tophit_and_contaminants2.py ${afanc_json} ${kraken_json} ${refseq} ${params.species} ${params.unmix_myco} ${resources} ${prev_species_json} ${params.permissive} ${pass}
    

    contam_to_remove=\$(jq -r '.summary_questions.are_there_contaminants' ${species_in_sample})
    acceptable_species=\$(jq -r '.summary_questions.is_the_top_species_appropriate' ${species_in_sample})
    top_hit=\$(jq -r '.top_hit.name' ${species_in_sample})

    if [ \$contam_to_remove == 'yes' ]; then
        if [ "${params.permissive}" == "no" ]; then
            printf "${sample_name}"
            echo '{"error":"sample remains contaminated, even after attempting to resolve this"}' | jq '.' > ${error_log} && jq -s ".[0] * .[1] * .[2]" ${software_json} ${error_log} ${species_in_sample} > ${report_json}
        else
            if [ "${pass}" == 2 ]; then
                 printf "NOW_ALIGN_TO_REF_${sample_name}"
            else
                 printf "${sample_name}"
            fi
            echo '{"warning":"sample remains contaminated, even after attempting to resolve this"}' | jq '.' > ${error_log} && jq -s ".[0] * .[1] * .[2]" ${software_json} ${error_log} ${species_in_sample} > ${report_json}
        fi
    fi

    if [ \$contam_to_remove == 'no' ] && [ \$acceptable_species == 'yes' ]; then 
        printf "NOW_ALIGN_TO_REF_${sample_name}"
    elif [ \$contam_to_remove == 'no' ] && [ \$acceptable_species == 'no' ]; then
        jq -n --arg key "\$top_hit" '{"error": ("top hit " + \$key + " does not have a reference genome. Sample will not proceed beyond preprocessing workflow.")}' > ${error_log} && \
        jq -s ".[0] * .[1] * .[2]" ${software_json} ${error_log} ${species_in_sample} > ${report_json}
        printf "DO_NOT_PROCEED_${sample_name}"
    fi
    """

    stub:
    error_log = "${sample_name}_err.json"

    """
    touch ${sample_name}_species_in_sample_pass_${pass}.json
    touch ${error_log}
    printf ${params.summary_doWeAlign}
    """
}

process count_pass {
    /**
    * @QCcheckpoint none
    */

    tag "${sample_name } pass ${pass}"
    label 'preprocessing'
    label 'low_cpu'
    label 'low_memory'

    input:
    val(sample_name) //just for tagging's sake
    val(pass)

    output:
    stdout emit: new_pass_number

    script:
    """
    echo \$(($pass + 1)) | tr -d [:space:]
    """
}
