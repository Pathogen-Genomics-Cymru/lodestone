// modules for the preprocessing workflow

process checkBamValidity {
    /**
    * @QCcheckpoint confirm that samtools validates bam
    */

    tag { bam_file.getBaseName() }
    label 'preprocessing'
    label 'low_memory'
    label 'low_cpu'
    
    publishDir "${params.output_dir}/${bam_file.getBaseName()}", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple path(bam_file), path(software_json)

    output:
    tuple path(bam_file), stdout, path(software_json), emit: checkValidity_bam
    path "${bam_file.getBaseName()}_err.json", emit: checkValidityBam_log optional true
    path "${bam_file.getBaseName()}_report.json", emit: checkValidityBam_report optional true

    script:
    error_log = "${bam_file.getBaseName()}_err.json"
    report_json = "${bam_file.getBaseName()}_report.json"

    """
    is_ok=\$(samtools quickcheck $bam_file && echo 'OK' || echo 'FAIL' )

    if [ \$is_ok == 'OK' ]; then printf \$is_ok; else echo '{"error":"bam did not pass samtools quickcheck"}' | jq '.' > ${error_log} && jq -s ".[0] * .[1]" ${software_json} ${error_log} > ${report_json}; fi
    """

    stub:
    error_log = "${bam_file.getBaseName()}_err.json"

    """
    touch ${error_log}
    printf ${params.checkBamValidity_isok}
    """
}

process checkFqValidity {
    /**
    * @QCcheckpoint confirm that fqtools validates both fastqs
    */

    tag { sample_name }
    label 'preprocessing'
    label 'low_memory'
    label 'low_cpu'

    errorStrategy 'ignore'

    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), path(software_json)

    output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, path(software_json), emit: checkValidity_fqs
    path "${sample_name}_err.json", emit: checkValidity_log optional true
    path "${sample_name}_report.json", emit: checkValidity_report optional true

    script:
    error_log = "${sample_name}_err.json"
    report_json = "${sample_name}_report.json"

    """
    is_ok=\$(fqtools validate $fq1 $fq2)

    if [ \$is_ok == 'OK' ]; then printf 'OK'; else echo '{"error":"sample did not pass fqtools validation check"}' | jq '.' > ${error_log} && jq -s ".[0] * .[1]" ${software_json} ${error_log} > ${report_json}; fi
    """

    stub:
    error_log  = "${sample_name}_err.json"

    """
    printf ${params.checkFqValidity_isok}
    touch ${error_log}
    """
}


process bam2fastq {
    /**
    * @QCcheckpoint none
    */

    tag { bam_file.getBaseName() }
    label 'preprocessing'
    label 'low_memory'
    label 'normal_cpu'

    input:
    tuple path(bam_file), val(is_ok), path(software_json)

    when:
    is_ok == 'OK'

    output:
    tuple val("${bam_file.getBaseName()}"), path("${bam_file.getBaseName()}_1.fq.gz"), path("${bam_file.getBaseName()}_2.fq.gz"), stdout, path(software_json), emit: bam2fastq_fqs

    script:
    """
    samtools sort -n $bam_file -o ${bam_file.getBaseName()}.sorted.bam

    bedtools bamtofastq -i ${bam_file.getBaseName()}.sorted.bam -fq ${bam_file.getBaseName()}_1.fq -fq2 ${bam_file.getBaseName()}_2.fq

    rm ${bam_file.getBaseName()}.sorted.bam

    gzip ${bam_file.getBaseName()}_1.fq || true
    gzip ${bam_file.getBaseName()}_2.fq || true

    printf 'OK'
    """

    stub:
    """
    touch ${bam_file.getBaseName()}_1.fq.gz
    touch ${bam_file.getBaseName()}_2.fq.gz
    printf 'OK'
    """
}

process countReads {
    /**
    * @QCcheckpoint fail sample if there are < 100k raw reads
    */

    tag { sample_name }
    label 'preprocessing'
    label 'low_memory'
    label 'low_cpu'

    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(is_ok), path(software_json)

    when:
    is_ok == 'OK'

    output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, path(software_json), emit: countReads_fqs
    path "${sample_name}_err.json", emit: countReads_log optional true
    path "${sample_name}_report.json", emit: countReads_report optional true

    script:
    error_log = "${sample_name}_err.json"
    report_json = "${sample_name}_report.json"

    """
    num_reads=\$(fqtools count $fq1 $fq2)

    if (( \$num_reads > 100000 )); then printf "${sample_name}"; else jq -n --arg key "\$num_reads" '{"error": ("sample did not have > 100k pairs of raw reads it only contained " + \$key)}' > ${error_log} && printf "fail" && jq -s ".[0] * .[1]" ${software_json} ${error_log} > ${report_json}; fi
    """

    stub:
    error_log = "${sample_name}_err.json"

    """
    printf ${params.countReads_runfastp}
    touch ${error_log}
    """
}

process fastp {
    /**
    * @QCcheckpoint confirm that there > 100k reads after cleaning with fastp
    */

    tag { sample_name }
    label 'preprocessing'
    label 'low_memory'
    label 'low_cpu'

    publishDir "${params.output_dir}/$sample_name/raw_read_QC_reports", mode: 'copy', pattern: '*_fastp.json'
    publishDir "${params.output_dir}/$sample_name/output_reads", mode: 'copy', pattern: '*.fq.gz' // may be overwritten if unmixing needed
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(run_fastp), path(software_json)

    when:
    run_fastp =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_cleaned_1.fq.gz"), path("${sample_name}_cleaned_2.fq.gz"), stdout, path(software_json), emit: fastp_fqs
    path("${sample_name}_fastp.json", emit: fastp_json)
    path "${sample_name}_err.json", emit: fastp_log optional true
    path "${sample_name}_report.json", emit: fastp_report optional true

    script:
    clean_fq1  = "${sample_name}_cleaned_1.fq.gz"
    clean_fq2  = "${sample_name}_cleaned_2.fq.gz"
    fastp_json = "${sample_name}_fastp.json"
    fastp_html = "${sample_name}_fastp.html"
    error_log  = "${sample_name}_err.json"
    report_json = "${sample_name}_report.json"

    """
    fastp -i $fq1 -I $fq2 -o ${clean_fq1} -O ${clean_fq2} -j ${fastp_json} -h ${fastp_html} --length_required 50 --average_qual 10 --low_complexity_filter --correction --cut_right --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20

    rm -rf ${fastp_html}

    num_reads=\$(fqtools count ${clean_fq1} ${clean_fq2})

    if (( \$num_reads > 100000 )); then printf "${sample_name}"; else jq -n --arg key "\$num_reads" '{"error": ("after fastp sample did not have > 100k pairs of raw reads it only contained " + \$key)}' > ${error_log} && printf "fail" && jq -s ".[0] * .[1]" ${software_json} ${error_log} > ${report_json}; fi
    """

    stub:
    clean_fq1  = "${sample_name}_cleaned_1.fq.gz"
    clean_fq2  = "${sample_name}_cleaned_2.fq.gz"
    fastp_json = "${sample_name}_fastp.json"
    fastp_html = "${sample_name}_fastp.html"
    error_log  = "${sample_name}_err.json"

    """
    printf ${params.fastp_enoughreads}
    touch ${error_log}
    touch ${clean_fq1}
    touch ${clean_fq2}
    touch ${fastp_json}
    touch ${fastp_html}
    """
}

process fastQC {
    /**
    * @QCcheckpoint none
    */

    tag { sample_name }
    label 'preprocessing'
    label 'low_memory'
    label 'low_cpu'

    publishDir "${params.output_dir}/$sample_name/raw_read_QC_reports", mode: 'copy'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_reads), path(software_json)

    output:
    path("*", emit: fastQC_all)

    script:
    """
    cat $fq1 $fq2 > ${sample_name}.fq.gz
    fastqc ${sample_name}.fq.gz
    rm ${sample_name}.fq.gz
    """

    stub:
    """
    touch ${sample_name}_fastqc.html
    touch ${sample_name}_fastqc.zip
    """
}

process kraken2 {
    /**
    * @QCcheckpoint if Kraken's top family classification is NOT Mycobacteriaceae, sample will not proceed further than afanc
    */

    tag { sample_name }
    label 'preprocessing'
    label 'normal_cpu'
    label 'high_memory'

    publishDir "${params.output_dir}/$sample_name/output_reads", mode: 'copy', pattern: '*.fq.gz', overwrite: 'true'
    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP", mode: 'copy', pattern: '*_kraken_report.*'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_reads), path(software_json)
    path(database)

    when:
    enough_reads =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_kraken_report.txt"), path("${sample_name}_kraken_report.json"), emit: kraken2_json
    tuple val(sample_name), path("${sample_name}_cleaned_1.fq.gz"), path("${sample_name}_cleaned_2.fq.gz"), stdout, path(software_json), emit: kraken2_fqs
    path "${sample_name}_err.json", emit: kraken2_log optional true
    path "${sample_name}_report.json", emit: kraken2_report optional true

    script:
    kraken2_report = "${sample_name}_kraken_report.txt"
    kraken2_json = "${sample_name}_kraken_report.json"
    kraken2_read_classification = "${sample_name}_read_classifications.txt"
    nonBac_depleted_reads_1 = "${sample_name}_cleaned_1.fq"
    nonBac_depleted_reads_2 = "${sample_name}_cleaned_2.fq"
    error_log = "${sample_name}_err.json"
    report_json = "${sample_name}_report.json"

    """
    kraken2 --threads ${task.cpus} --db . --output ${kraken2_read_classification} --report ${kraken2_report} --paired $fq1 $fq2
    
    parse_kraken_report2.py ${kraken2_report} ${kraken2_json} ${params.percent_threshold} ${params.n_reads_threshold}

    extract_kraken_reads.py -k ${kraken2_read_classification} -r ${kraken2_report} -s $fq1 -s2 $fq2 -o ${nonBac_depleted_reads_1} -o2 ${nonBac_depleted_reads_2} --taxid 2 --include-children --fastq-output >/dev/null

    gzip -f ${nonBac_depleted_reads_1}
    gzip -f ${nonBac_depleted_reads_2}

    rm -rf ${sample_name}_read_classifications.txt

    run_afanc=\$(jq '.afanc' ${kraken2_json})

    if [ \$run_afanc == '\"true\"' ]; then printf "${sample_name}"; else echo '{"error":"Kraken's top family hit either wasn't Mycobacteriaceae, or there were < 100k Mycobacteriaceae reads. Sample will not proceed further than afanc."}' | jq '.' > ${error_log} && printf "no" && jq -s ".[0] * .[1]" ${software_json} ${error_log} > ${report_json}; fi
    """

    stub:
    kraken2_report = "${sample_name}_kraken_report.txt"
    kraken2_json = "${sample_name}_kraken_report.json"
    kraken2_read_classification = "${sample_name}_read_classifications.txt"
    nonBac_depleted_reads_1 = "${sample_name}_cleaned_1.fq.gz"
    nonBac_depleted_reads_2 = "${sample_name}_cleaned_2.fq.gz"
    error_log = "${sample_name}_err.json"

    """
    printf ${params.kraken2_runmykrobe}
    touch ${kraken2_report}
    touch ${kraken2_json}
    touch ${kraken2_read_classification}
    touch ${nonBac_depleted_reads_1}
    touch ${nonBac_depleted_reads_2}
    touch ${error_log}
    """
}

process afanc {
    /**
    * @QCcheckpoint none
    */

    tag { sample_name }
    label 'preprocessing'
    label 'normal_cpu'
    label 'high_memory'
    label 'retry_afanc'

    errorStrategy 'ignore'

    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP", mode: 'copy', pattern: '*_afanc*.json'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(run_afanc), path(software_json), path(kraken_report), path(kraken_json)
    path(afanc_myco_db)
    val(resource_dir)
    path(refseq_path)

    output:
    tuple val(sample_name), path("${sample_name}_afanc_report.json"), stdout, emit: afanc_json
    path "${sample_name}_err.json", emit: afanc_log optional true
    path "${sample_name}_report.json", emit: afanc_report optional true
    path "${sample_name}_afanc_original.json", emit: afanc_original optional true

    script:
    afanc_report = "${sample_name}_afanc_report.json"
    error_log = "${sample_name}_err.json"
    report_json = "${sample_name}_report.json"

    """
    if [[ ${run_afanc} =~ /${sample_name}/ ]]
    then
	afanc screen ${afanc_myco_db} ${fq1} ${fq2} -p 5.0 -n 1000 -o ${sample_name} -t ${task.cpus} -v ${afanc_myco_db}/lineage_profiles/TB_variants.tsv
        cp ${sample_name}/${sample_name}.json ${sample_name}_afanc_original.json
	reformat_afanc_json.py ${sample_name}/${sample_name}.json
	printf ${sample_name}
    else
	afanc screen ${afanc_myco_db} ${fq1} ${fq2} -p 2.0 -n 500 -o ${sample_name} -t ${task.cpus} -v ${afanc_myco_db}/lineage_profiles/TB_variants.tsv
        cp ${sample_name}/${sample_name}.json ${sample_name}_afanc_original.json
	reformat_afanc_json.py ${sample_name}/${sample_name}.json

	identify_tophit_and_contaminants2.py ${afanc_report} ${kraken_json} $refseq_path ${params.species} ${params.unmix_myco} $resource_dir null

	echo '{"error":"Kraken's top family hit either wasn't Mycobacteriaceae, or there were < 100k Mycobacteriaceae reads. Sample will not proceed further than afanc."}' | jq '.' > ${error_log} && printf "no" && jq -s ".[0] * .[1] * .[2]" ${software_json} ${error_log} ${sample_name}_species_in_sample.json > ${report_json}

    fi

    """

    stub:
    afanc_report = "${sample_name}_afanc_report.json"

    """
    touch ${afanc_report}
    printf ${sample_name}
    """
}


process mykrobe {
    /**
    * @QCcheckpoint none
    */

    tag { sample_name }
    label 'preprocessing'
    label 'normal_cpu'
    label 'medium_memory'

    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP", mode: 'copy', pattern: '*_mykrobe_report.json'
    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP", mode: 'copy', pattern: '*_mykrobe_report.csv'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(run_mykrobe), path(software_json)

    when:
    run_mykrobe =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_mykrobe_report.json"), stdout, emit: mykrobe_report
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: mykrobe_fqs

    script:
    mykrobe_report = "${sample_name}_mykrobe_report"

    """
    mykrobe predict --sample ${sample_name} --species tb --threads ${task.cpus} --format json_and_csv --output ${mykrobe_report} -1 $fq1 $fq2
    printf ${sample_name}
    """

    stub:
    mykrobe_report = "${sample_name}_mykrobe_report.json"

    """
    touch ${mykrobe_report}
    printf ${sample_name}
    """
}

process bowtie2 {
    /**
    * @QCcheckpoint none
    */

    tag { sample_name }
    label 'preprocessing'
    label 'normal_cpu'
    label 'medium_memory'

    publishDir "${params.output_dir}/$sample_name/output_reads", mode: 'copy', pattern: '*.fq.gz', overwrite: 'true'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_myco_reads), path(software_json)
    path(index)

    when:
    enough_myco_reads =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_cleaned_1.fq.gz"), path("${sample_name}_cleaned_2.fq.gz"), path(software_json), emit: bowtie2_fqs

    script:
    bam = "${sample_name}.bam"
    humanfree_fq1 = "${sample_name}_cleaned_1.fq"
    humanfree_fq2 = "${sample_name}_cleaned_2.fq"

    """
    bowtie2 --very-sensitive -p ${task.cpus} -x ./${params.bowtie_index_name} -1 $fq1 -2 $fq2 | samtools view -f 4 -Shb - > ${bam}
    samtools fastq -1 ${humanfree_fq1} -2 ${humanfree_fq2} -s singleton.fq ${bam}

    rm -rf ${bam}
    rm -rf singleton.fq

    gzip -f ${humanfree_fq1}
    gzip -f ${humanfree_fq2}
    """

    stub:
    humanfree_fq1 = "${sample_name}_cleaned_1.fq"
    humanfree_fq2 = "${sample_name}_cleaned_2.fq"

    """
    touch ${humanfree_fq1}.gz
    touch ${humanfree_fq2}.gz
    """
}

process identifyBacterialContaminants {
    /**
    * @QCcheckpoint if urllist.txt is empty, there are no contaminant genomes to download, so skip next process
    */

    tag { sample_name }
    label 'preprocessing'
    label 'normal_cpu'
    label 'medium_memory'

    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP", mode: 'copy', pattern: '*_species_in_samp*.json'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), path(software_json), path(afanc_json), val(enough_myco_reads), path(kraken_report), path(kraken_json)
    val(resources)
    path(refseq)

    when:
    enough_myco_reads =~ /${sample_name}/

    output:// enable dsl2
nextflow.enable.dsl = 2
    tuple val(sample_name), path("${sample_name}_urllist.txt"), stdout, path(software_json), path("${sample_name}_species_in_sample_previous.json"), emit: contam_list optional true
    tuple val(sample_name), path("${sample_name}_species_in_sample_previous.json"), stdout, path(software_json), emit: prev_sample_json optional true
    tuple val(sample_name), path("${sample_name}_species_in_sample.json"), stdout, emit: sample_json
    tuple val(sample_name), path("${sample_name}_nocontam_1.fq.gz"), path("${sample_name}_nocontam_2.fq.gz"), path(software_json), path("${sample_name}_species_in_sample.json"), stdout, emit: nocontam_fqs optional true
    path "${sample_name}_err.json", emit: contam_log optional true
    path "${sample_name}_report.json", emit: contam_report optional true

    script:
    error_log = "${sample_name}_err.json"
    report_json = "${sample_name}_report.json"

    """
    identify_tophit_and_contaminants2.py ${afanc_json} ${kraken_json} ${refseq} ${params.species} ${params.unmix_myco} ${resources} null

    contam_to_remove=\$(jq -r '.summary_questions.are_there_contaminants' ${sample_name}_species_in_sample.json)
    acceptable_species=\$(jq -r '.summary_questions.is_the_top_species_appropriate' ${sample_name}_species_in_sample.json)
    top_hit=\$(jq -r '.top_hit.name' ${sample_name}_species_in_sample.json)

    if [ \$contam_to_remove == 'yes' ]; then cp ${sample_name}_species_in_sample.json ${sample_name}_species_in_sample_previous.json; fi

    if [ \$contam_to_remove == 'yes' ]; then printf "NOW_DECONTAMINATE_${sample_name}"; elif [ \$contam_to_remove == 'no' ] && [ \$acceptable_species == 'yes' ]; then printf "NOW_ALIGN_TO_REF_${sample_name}" && mv $fq1 ${sample_name}_nocontam_1.fq.gz && mv $fq2 ${sample_name}_nocontam_2.fq.gz; elif [ \$contam_to_remove == 'no' ] && [ \$acceptable_species == 'no' ]; then jq -n --arg key "\$top_hit" '{"error": ("top hit " + \$key + " does not have a reference genome. Sample will not proceed beyond preprocessing workflow.")}' > ${error_log} && jq -s ".[0] * .[1] * .[2]" ${software_json} ${error_log} ${sample_name}_species_in_sample.json > ${report_json}; fi
    """

    stub:
    error_log = "${sample_name}_err.json"

    """
    touch ${sample_name}_species_in_sample.json
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

    tag { sample_name }
    label 'preprocessing'
    label 'low_cpu'
    label 'medium_memory'

    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(contam_list), val(run_decontaminator), path(software_json), path(prev_species_json)

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

    tag { sample_name }
    label 'preprocessing'
    label 'normal_cpu'
    label 'high_memory'

    publishDir "${params.output_dir}/$sample_name/output_reads", mode: 'copy', pattern: '*.fq.gz', overwrite: 'true'

    input:
    tuple val(sample_name), path(fq1), path(fq2), path(software_json), path(contam_fa), val(does_fa_pass)

    when:
    does_fa_pass =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_cleaned_1.fq.gz"), path("${sample_name}_cleaned_2.fq.gz"), path(software_json), emit: reClassification_fqs

    script:
    bam = "${sample_name}.bam"
    decontam_fq1 = "${sample_name}_cleaned_1.fq"
    decontam_fq2 = "${sample_name}_cleaned_2.fq"

    """
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

    tag { sample_name }
    label 'preprocessing'
    label 'normal_cpu'
    label 'high_memory'

    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP_and_postContamRemoval", mode: 'copy', pattern: '*_kraken*'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(fq1), path(fq2), path(software_json)
    path(database)

    output:
    tuple val(sample_name), path("${sample_name}_kraken_report.txt"), path("${sample_name}_kraken_report.json"), emit: reKraken_report

    script:
    kraken2_report = "${sample_name}_kraken_report.txt"
    kraken2_json = "${sample_name}_kraken_report.json"
    kraken2_read_classification = "${sample_name}_read_classifications.txt"

    """
    kraken2 --threads ${task.cpus} --db . --output ${kraken2_read_classification} --report ${kraken2_report} --paired $fq1 $fq2

    parse_kraken_report2.py ${kraken2_report} ${kraken2_json} ${params.percent_threshold} ${params.n_reads_threshold}
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

    tag { sample_name }
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

    output:
    tuple val(sample_name), path("${sample_name}_afanc_report.json"), emit: reAfanc_report
    path "${sample_name}_afanc_original.json", emit: reAfanc_original
    // tuple val(sample_name), path(fq1), path(fq2), stdout, emit: afanc_fqs

    script:
    afanc_report = "${sample_name}_afanc_report.json"

    """
    afanc screen ${afanc_myco_db} ${fq1} ${fq2} -p 5.0 -n 1000 -o ${sample_name} -t ${task.cpus} -v ${afanc_myco_db}/lineage_profiles/TB_variants.tsv

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

    tag { sample_name }
    label 'preprocessing'
    label 'normal_cpu'
    label 'low_memory'

    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP_and_postContamRemoval", mode: 'copy', pattern: '*_mykrobe_report.json'
    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP_and_postContamRemoval", mode: 'copy', pattern: '*_mykrobe_report.csv'

    input:
    tuple val(sample_name), path(fq1), path(fq2), path(software_json)

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

    tag { sample_name }
    label 'preprocessing'
    label 'low_cpu'
    label 'medium_memory'

    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP_and_postContamRemoval", mode: 'copy', pattern: '*_species_in_sample.json'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*{_err.json,_report.json}'

    input:
    tuple val(sample_name), path(afanc_json), path(kraken_report), path(kraken_json), path(prev_species_json), val(decontam), path(software_json)
    val(resources)
    path(refseq)

    output:
    tuple val(sample_name), path("${sample_name}_species_in_sample.json"), stdout, emit: summary_json
    path "${sample_name}_err.json", emit: summary_log optional true
    path "${sample_name}_report.json", emit: summary_report optional true

    script:
    error_log = "${sample_name}_err.json"
    report_json = "${sample_name}_report.json"

    """
    identify_tophit_and_contaminants2.py ${afanc_json} ${kraken_json} ${refseq} ${params.species} ${params.unmix_myco} ${resources} ${prev_species_json}
    

    contam_to_remove=\$(jq -r '.summary_questions.are_there_contaminants' ${sample_name}_species_in_sample.json)
    acceptable_species=\$(jq -r '.summary_questions.is_the_top_species_appropriate' ${sample_name}_species_in_sample.json)
    top_hit=\$(jq -r '.top_hit.name' ${sample_name}_species_in_sample.json)

    if [ \$contam_to_remove == 'yes' ]; then echo '{"error":"sample remains contaminated, even after attempting to resolve this"}' | jq '.' > ${error_log} && jq -s ".[0] * .[1] * .[2]" ${software_json} ${error_log} ${sample_name}_species_in_sample.json > ${report_json}; fi

    if [ \$contam_to_remove == 'no' ] && [ \$acceptable_species == 'yes' ]; then printf "NOW_ALIGN_TO_REF_${sample_name}"; elif [ \$contam_to_remove == 'no' ] && [ \$acceptable_species == 'no' ]; then jq -n --arg key "\$top_hit" '{"error": ("top hit " + \$key + " does not have a reference genome. Sample will not proceed beyond preprocessing workflow.")}' > ${error_log} && jq -s ".[0] * .[1] * .[2]" ${software_json} ${error_log} ${sample_name}_species_in_sample.json > ${report_json}; fi
    """

    stub:
    error_log = "${sample_name}_err.json"

    """
    touch ${sample_name}_species_in_sample.json
    touch ${error_log}
    printf ${params.summary_doWeAlign}
    """
}
