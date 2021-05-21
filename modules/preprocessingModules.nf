// modules for the preprocessing workflow

process checkBamValidity {
    /**
    * @QCcheckpoint confirm that samtools validates bam
    */

    tag { bam_file.getBaseName() }
  
    publishDir "${params.output_dir}/${bam_file.getBaseName()}", mode: 'copy', pattern: '*.log'

    memory '5 GB'

    input:
    path(bam_file)
   
    output:
    tuple path(bam_file), stdout, emit: checkValidity_bam

    script:
    error_log = "${bam_file.getBaseName()}.log"
			
    """
    is_ok=\$(samtools quickcheck $bam_file && echo 'OK' || echo 'FAIL' )

    if [ \$is_ok == 'OK' ]; then printf "" >> ${error_log}; else echo "error: bam did not pass samtools quickcheck" >> ${error_log}; fi
    printf \$is_ok
    """
}

process checkFqValidity {
    /**
    * @QCcheckpoint confirm that fqtools validates both fastqs
    */
   
    tag { sample_name }

    publishDir "${params.output_dir}/$sample_name", mode: 'copy', pattern: '*.log'

    memory '5 GB'

    input:
    tuple val(sample_name), path(fq1), path(fq2)
	
    output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: checkValidity_fqs
    path("${sample_name}.log", emit: checkValidity_log)
		
    script:
    error_log  = "${sample_name}.log"
	
    """
    is_ok=\$(fqtools validate $fq1 $fq2)

    if [ \$is_ok == 'OK' ]; then printf "" >> ${error_log}; else echo "error: sample did not pass fqtools validation check" >> ${error_log}; fi
    printf \$is_ok
    """
}


process bam2fastq {
    /**
    * @QCcheckpoint none
    */

    tag { bam_file.getBaseName() }

    memory '5 GB'
    
    when:
    is_ok == 'OK'    

    input:
    path(bam_file)
    
    output:
    tuple val("${bam_file.getBaseName()}"), path("${bam_file.getBaseName()}_1.fq.gz"), path("${bam_file.getBaseName()}_2.fq.gz"), emit: bam2fastq_fqs

    script:
    """
    samtools sort -n $bam_file -o ${bam_file.getBaseName()}.sorted.bam

    bedtools bamtofastq -i ${bam_file.getBaseName()}.sorted.bam -fq ${bam_file.getBaseName()}_1.fq -fq2 ${bam_file.getBaseName()}_2.fq

    rm ${bam_file.getBaseName()}.sorted.bam

    gzip ${bam_file.getBaseName()}_1.fq || true
    gzip ${bam_file.getBaseName()}_2.fq || true

    printf 'OK'
    """
}

process countReads {
    /**
    * @QCcheckpoint fail sample if there are < 100k raw reads
    */

    tag { sample_name }

    publishDir "${params.output_dir}/$sample_name", mode: 'copy', pattern: '*.log'

    memory '5 GB'

    when:
    is_ok == 'OK'
  
    input:
    tuple val(sample_name), path(fq1), path(fq2), val(is_ok)

    output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: countReads_fqs
    path("${sample_name}.log", emit: countReads_log)

    script:
    error_log = "${sample_name}.log"
    """
    num_reads=\$(fqtools count $fq1 $fq2)

    if (( \$num_reads > 100000 )); then printf "" >> ${error_log} && printf "pass"; else echo "error: sample did not have > 100k pairs of raw reads (it only contained \$num_reads)" >> ${error_log} && printf "fail"; fi
    """
}
	
process fastp {
    /**
    * @QCcheckpoint confirm that there > 100k reads after cleaning with fastp
    */
     
    tag { sample_name }
 
    publishDir "${params.output_dir}/$sample_name/raw_read_QC_reports", mode: 'copy', pattern: '*.json'
    publishDir "${params.output_dir}/$sample_name/output_reads", mode: 'copy', pattern: '*.fq.gz' // may be overwritten if unmixing needed
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', pattern: '*.log'

    memory '5 GB'

    when:
    run_fastp == 'pass'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(run_fastp)
		
    output:
    tuple val(sample_name), path("${sample_name}_cleaned_1.fq.gz"), path("${sample_name}_cleaned_2.fq.gz"), stdout, emit: fastp_fqs
    path("${sample_name}_fastp.json", emit: fastp_json)
    path("${sample_name}.log", emit: fastp_log)
   
    script:
    clean_fq1  = "${sample_name}_cleaned_1.fq.gz"
    clean_fq2  = "${sample_name}_cleaned_2.fq.gz"
    fastp_json = "${sample_name}_fastp.json"
    fastp_html = "${sample_name}_fastp.html"
    error_log  = "${sample_name}.log"
	
    """
    fastp -i $fq1 -I $fq2 -o ${clean_fq1} -O ${clean_fq2} -j ${fastp_json} -h ${fastp_html} --length_required 50 --average_qual 10 --low_complexity_filter --correction --cut_right --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20
    
    rm -rf ${fastp_html}

    num_reads=\$(jq '.summary.after_filtering.total_reads' ${fastp_json} | awk '{sum+=\$0} END{print sum}')

    if (( \$num_reads > 100000 )); then printf "" >> ${error_log} && printf "pass"; else echo "error: after fastp, sample did not have > 100k reads (it only contained \$num_reads)" >> ${error_log} && printf "fail"; fi
    """
}

process fastQC {
    /**
    * @QCcheckpoint none
    */
	
    tag { sample_name }

    publishDir "${params.output_dir}/$sample_name/raw_read_QC_reports", mode: 'copy'

    memory '5 GB'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_reads)
	
    output:
    path("*", emit: fastQC_all)
	
    script:
    """
    cat $fq1 $fq2 > ${sample_name}.fq.gz
    fastqc ${sample_name}.fq.gz
    rm ${sample_name}.fq.gz
    """
}

process kraken2 {
    /**
    * @QCcheckpoint only pass to Mykrobe if Kraken's top family classification is Mycobacteriaceae
    */

    tag { sample_name }

    publishDir "${params.output_dir}/$sample_name/speciation_reports_cleanedReads", mode: 'copy', pattern: '*_kraken_report.*'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', pattern: '*.log'
       
    cpus 8

    memory '10 GB'

    when:
    enough_reads == 'pass'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_reads)
    path(database)
		
    output:
    tuple path("${sample_name}_kraken_report.txt"), path("${sample_name}_kraken_report.json"), emit: kraken2_report
    tuple val(sample_name), path("${sample_name}_clean_exclNonBacteria_1.fq.gz"), path("${sample_name}_clean_exclNonBacteria_2.fq.gz"), stdout, emit: kraken2_fqs
    path("${sample_name}.log", emit: kraken2_log)
			
    script:
    kraken2_report = "${sample_name}_kraken_report.txt"
    kraken2_json = "${sample_name}_kraken_report.json"
    kraken2_read_classification = "${sample_name}_read_classifications.txt"
    nonBac_depleted_reads_1 = "${sample_name}_clean_exclNonBacteria_1.fq"
    nonBac_depleted_reads_2 = "${sample_name}_clean_exclNonBacteria_2.fq"
    error_log = "${sample_name}.log"
	
    """
    kraken2 --threads ${task.cpus} --db . --output ${kraken2_read_classification} --report ${kraken2_report} --paired $fq1 $fq2

    perl ${baseDir}/bin/parse_kraken_report2.pl ${kraken2_report} ${kraken2_json} 0.5 5000

    ${baseDir}/bin/extract_kraken_reads.py -k ${kraken2_read_classification} -r ${kraken2_report} -s $fq1 -s2 $fq2 -o ${nonBac_depleted_reads_1} -o2 ${nonBac_depleted_reads_2} --taxid 2 --include-children --fastq-output >/dev/null

    gzip ${nonBac_depleted_reads_1}
    gzip ${nonBac_depleted_reads_2}

    rm -rf ${sample_name}_read_classifications.txt

    run_mykrobe=\$(jq '.Mykrobe' ${kraken2_json})

    if [ \$run_mykrobe == '\"true\"' ]; then printf "" >> ${error_log} && printf "yes"; else echo "error: Kraken's top family hit either wasn't Mycobacteriaceae, or represented < 5000 reads in absolute terms or < 0.5% of the total" >> ${error_log} && printf "no"; fi
    """
}

process mykrobe {
    /**
    * @QCcheckpoint none
    */

    tag { sample_name }

    publishDir "${params.output_dir}/$sample_name/speciation_reports_cleanedReads", mode: 'copy', pattern: '*_mykrobe_report.json'

    cpus 8

    memory '5 GB'

    when:
    run_mykrobe == 'yes'
	
    input:
    tuple val(sample_name), path(fq1), path(fq2), val(run_mykrobe)
		
    output:
    tuple val(sample_name), path("${sample_name}_mykrobe_report.json"), stdout, emit: mykrobe_report

    script:
    mykrobe_report = "${sample_name}_mykrobe_report.json"
	
    """
    mykrobe predict ${sample_name} tb --threads ${task.cpus} --format json --output ${mykrobe_report} -1 $fq1 $fq2
    printf 'yes'
    """
}

process bowtie2 {
    /**
    * @QCcheckpoint none
    */

    tag { sample_name }

    publishDir "${params.output_dir}/$sample_name/output_reads", mode: 'copy', pattern: '*.fq.gz', overwrite: 'true'

    cpus 8

    memory '5 GB'

    when:
    enough_myco_reads == 'yes'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_myco_reads)
    path(index)

    output:
    tuple val(sample_name), path("${sample_name}_cleaned_1.fq.gz"), path("${sample_name}_cleaned_2.fq.gz"), emit: bowtie2_fqs

    script:
    bam = "${sample_name}.bam"
    humanfree_fq1 = "${sample_name}_cleaned_1.fq"
    humanfree_fq2 = "${sample_name}_cleaned_2.fq"
	
    """
    bowtie2 --very-sensitive -p ${task.cpus} -x ${index}/${params.bowtie_index_name} -1 $fq1 -2 $fq2 | samtools view -f 4 -Shb - > ${bam}
    samtools fastq -1 ${humanfree_fq1} -2 ${humanfree_fq2} -s singleton.fq ${bam}

    rm -rf ${bam}
    rm -rf singleton.fq

    gzip -f ${humanfree_fq1}
    gzip -f ${humanfree_fq2}
    """	
}

process identifyBacterialContaminants {
    /**
    * @QCcheckpoint if urllist.txt is empty, there are no contaminant genomes to download, so skip next process
    */    

    tag { sample_name }

    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP", mode: 'copy', pattern: '*_species_in_sample.json'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', pattern: '*.log' 

    when:
    enough_myco_reads == 'yes'

    input:
    tuple val(sample_name), path(mykrobe_json), val(enough_myco_reads)
    tuple path(kraken_report), path(kraken_json)
		
    output:
    tuple val(sample_name), path("${sample_name}_urllist.txt"), stdout, emit: contam_list
    path("${sample_name}_species_in_sample_previous.json"), emit: prev_sample_json
    path("${sample_name}.log", emit: contam_log)
		
    script:
    error_log = "${sample_name}.log"	

    """
    perl ${baseDir}/bin/identify_tophit_and_contaminants2.pl ${mykrobe_json} ${kraken_json} ${baseDir}/resources/assembly_summary_refseq.txt ${params.species} ${params.unmix_myco} ${baseDir}/resources null
    
    cp ${sample_name}_species_in_sample.json ${sample_name}_species_in_sample_previous.json

    contam_to_remove=\$(jq -r '.summary_questions.are_there_contaminants' ${sample_name}_species_in_sample.json)
    acceptable_species=\$(jq -r '.summary_questions.is_the_top_species_appropriate' ${sample_name}_species_in_sample.json)
    top_hit=\$(jq -r '.top_hit.name' ${sample_name}_species_in_sample.json)

    if [ \$contam_to_remove == 'yes' ]; then printf "NOW_DECONTAMINATE_${sample_name}" && printf "" >> ${error_log}; elif [ \$contam_to_remove == 'no' ] && [ \$acceptable_species == 'yes' ]; then printf "NOW_ALIGN_TO_REF_${sample_name}" && printf "" >> ${error_log}; elif [ \$contam_to_remove == 'no' ] && [ \$acceptable_species == 'no' ]; then echo "top hit (\$top_hit) is not one of the 10 accepted mycobacteria" >> ${error_log}; fi
    """
}

process downloadContamGenomes {
    /**
    * @QCcheckpoint confirm that we could download every genome in the list of contaminants
    */
    
    tag { sample_name }

    publishDir "${params.output_dir}/$sample_name", mode: 'copy', pattern: '*.log'

    when:
    run_decontaminator =~ /NOW\_DECONTAMINATE\_${sample_name}/

    input:
    tuple val(sample_name), path(contam_list), val(run_decontaminator)
		
    output:
    tuple path("${sample_name}_contaminants.fa"), stdout, emit: contam_fa
    path("${sample_name}.log", emit: downcontam_log)
	
    script:
    contaminant_fa = "${sample_name}_contaminants.fa"
    error_log = "${sample_name}.log"
	
    """
    wget -i ${contam_list} --spider -nv -a linktestlog.txt 2>&1
    cat linktestlog.txt | awk '/listing\" \\[1\\]/{print \$4}' > confirmedurllist.txt
    wget -i confirmedurllist.txt

    gunzip *.gz
    cat *.fna > ${contaminant_fa}
    rm -rf *.fna

    num_urls_in=\$(cat $contam_list | wc -l)
    num_urls_out=\$(cat confirmedurllist.txt | wc -l)

    rm -rf linktestlog.txt confirmedurllist.txt

    if (( \$num_urls_in == \$num_urls_out )); then printf "" >> ${error_log} && printf "pass"; else echo "error: there were \$num_urls_in contaminant genomes but only \$num_urls_out could be downloaded" >> ${error_log} && printf "fail"; fi
    """
}

process mapToContamFa {
    /**
    * @QCcheckpoint none
    */

    tag { sample_name }

    publishDir "${params.output_dir}/$sample_name/output_reads", mode: 'copy', pattern: '*.fq.gz', overwrite: 'true'

    cpus 8

    memory '10 GB'

    when:
    does_fa_pass == 'pass'

    input:
    tuple val(sample_name), path(fq1), path(fq2)
    tuple path(contam_fa), val(does_fa_pass)
			
    output:
    tuple val(sample_name), path("${sample_name}_cleaned_1.fq.gz"), path("${sample_name}_cleaned_2.fq.gz"), emit: reClassification_fqs

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
}

process reKraken {
    /**
    * @QCcheckpoint none
    */

    tag { sample_name }

    publishDir "${params.output_dir}/$sample_name/speciation_reports_cleanedAndUnmixedReads", mode: 'copy', pattern: '*_kraken_report.*'

    cpus 8

    memory '10 GB'
    
    input:
    tuple val(sample_name), path(fq1), path(fq2)
    path(database)
		
    output:
    tuple path("${sample_name}_kraken_report.txt"), path("${sample_name}_kraken_report.json"), emit: reKraken_report

    script:
    kraken2_report = "${sample_name}_kraken_report.txt"
    kraken2_json = "${sample_name}_kraken_report.json"
    kraken2_read_classification = "${sample_name}_read_classifications.txt"
    
    """
    kraken2 --threads ${task.cpus} --db . --output ${kraken2_read_classification} --report ${kraken2_report} --paired $fq1 $fq2

    perl ${baseDir}/bin/parse_kraken_report2.pl ${kraken2_report} ${kraken2_json} 0.5 5000
    rm -rf ${sample_name}_read_classifications.txt
    """
}

process reMykrobe {
    /**
    * @QCcheckpoint none
    */
     
    tag { sample_name }

    publishDir "${params.output_dir}/$sample_name/speciation_reports_cleanedAndUnmixedReads", mode: 'copy', pattern: '*_mykrobe_report.json'
	
    cpus 8

    memory '5 GB'

    input:
    tuple val(sample_name), path(fq1), path(fq2)
		
    output:
    tuple val(sample_name), path("${sample_name}_mykrobe_report.json"), emit: reMykrobe_report

    script:
    mykrobe_report = "${sample_name}_mykrobe_report.json"
	
    """
    mykrobe predict ${sample_name} tb --threads ${task.cpus} --format json --output ${mykrobe_report} -1 $fq1 $fq2
    """
}

process summarise {
    /**
    * @QCcheckpoint none
    */    

    tag { sample_name }

    publishDir "${params.output_dir}/$sample_name/speciation_reports_cleanedAndUnmixedReads", mode: 'copy', pattern: '*.json'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', pattern: '*.log'

    input:
    tuple val(sample_name), path(mykrobe_json)
    tuple path(kraken_report), path(kraken_json)
    path(prev_species_json)
		
    output:
    path("${sample_name}_species_in_sample.json", emit: summary_json)
    path("${sample_name}.log", emit: summary_log)
	
    script:
    error_log = "${sample_name}.log"
	
    """
    perl ${baseDir}/bin/identify_tophit_and_contaminants2.pl ${mykrobe_json} ${kraken_json} ${baseDir}/resources/assembly_summary_refseq.txt ${params.species} ${params.unmix_myco} ${baseDir}/resources ${prev_species_json}
	
    contam_to_remove=\$(jq -r '.summary_questions.are_there_contaminants' ${sample_name}_species_in_sample.json)
    acceptable_species=\$(jq -r '.summary_questions.is_the_top_species_appropriate' ${sample_name}_species_in_sample.json)
    top_hit=\$(jq -r '.top_hit.name' ${sample_name}_species_in_sample.json)
	
    if [ \$contam_to_remove == 'yes' ]; then echo "error: sample remains contaminated, even after attempting to resolve this" >> ${error_log}; fi
	
    if [ \$contam_to_remove == 'no' ] && [ \$acceptable_species == 'yes' ]; then printf "NOW_ALIGN_TO_REF_${sample_name}" && printf "" >> ${error_log}; elif [ \$contam_to_remove == 'no' ] && [ \$acceptable_species == 'no' ]; then echo "error: top hit (\$top_hit) is not one of the 10 accepted mycobacteria" >> ${error_log}; fi
    """
}

