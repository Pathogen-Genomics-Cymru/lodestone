// modules for the preprocessing workflow

process checkBamValidity {
    /**
    * @QCcheckpoint confirm that samtools validates bam
    */

    tag { bam_file.getBaseName() }
    label 'preprocessing'
    label 'low_memory'

    publishDir "${params.output_dir}/${bam_file.getBaseName()}", mode: 'copy', overwrite: 'true', pattern: '*_err.json'

    input:
    path(bam_file)

    output:
    tuple path(bam_file), stdout, emit: checkValidity_bam

    script:
    error_log = "${bam_file.getBaseName()}_err.json"

    """
    is_ok=\$(samtools quickcheck $bam_file && echo 'OK' || echo 'FAIL' )

    if [ \$is_ok == 'OK' ]; then printf "" >> ${error_log}; else echo '{"error":"bam did not pass samtools quickcheck"}' | jq '.' >> ${error_log}; fi
    printf \$is_ok
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

    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*_err.json'

    input:
    tuple val(sample_name), path(fq1), path(fq2)

    output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: checkValidity_fqs
    path("${sample_name}_err.json", emit: checkValidity_log)

    script:
    error_log = "${sample_name}_err.json"

    """
    is_ok=\$(fqtools validate $fq1 $fq2)

    if [ \$is_ok == 'OK' ]; then printf 'OK' && printf "" >> ${error_log}; else echo '{"error":"sample did not pass fqtools validation check"}' | jq '.' >> ${error_log}; fi
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

    input:
    tuple path(bam_file), val(is_ok)

    when:
    is_ok == 'OK'

    output:
    tuple val("${bam_file.getBaseName()}"), path("${bam_file.getBaseName()}_1.fq.gz"), path("${bam_file.getBaseName()}_2.fq.gz"), stdout, emit: bam2fastq_fqs

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

    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*_err.json'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(is_ok)

    when:
    is_ok == 'OK'

    output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: countReads_fqs
    path("${sample_name}_err.json", emit: countReads_log)

    script:
    error_log = "${sample_name}_err.json"
    """
    num_reads=\$(fqtools count $fq1 $fq2)

    if (( \$num_reads > 100000 )); then printf "" >> ${error_log} && printf "${sample_name}"; else echo '{"error":"sample did not have > 100k pairs of raw reads (it only contained \$num_reads)"}' | jq '.' >> ${error_log} && printf "fail"; fi
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

    publishDir "${params.output_dir}/$sample_name/raw_read_QC_reports", mode: 'copy', pattern: '*_fastp.json'
    publishDir "${params.output_dir}/$sample_name/output_reads", mode: 'copy', pattern: '*.fq.gz' // may be overwritten if unmixing needed
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*_err.json'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(run_fastp)

    when:
    run_fastp =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_cleaned_1.fq.gz"), path("${sample_name}_cleaned_2.fq.gz"), stdout, emit: fastp_fqs
    path("${sample_name}_fastp.json", emit: fastp_json)
    path("${sample_name}_err.json", emit: fastp_log)

    script:
    clean_fq1  = "${sample_name}_cleaned_1.fq.gz"
    clean_fq2  = "${sample_name}_cleaned_2.fq.gz"
    fastp_json = "${sample_name}_fastp.json"
    fastp_html = "${sample_name}_fastp.html"
    error_log  = "${sample_name}_err.json"

    """
    fastp -i $fq1 -I $fq2 -o ${clean_fq1} -O ${clean_fq2} -j ${fastp_json} -h ${fastp_html} --length_required 50 --average_qual 10 --low_complexity_filter --correction --cut_right --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20

    rm -rf ${fastp_html}

    num_reads=\$(fqtools count $fq1 $fq2)

    if (( \$num_reads > 100000 )); then printf "" >> ${error_log} && printf "${sample_name}"; else echo '{"error":"after fastp, sample did not have > 100k pairs of reads (it only contained \$num_reads)"}' | jq '.' >> ${error_log} && printf "fail"; fi
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

    publishDir "${params.output_dir}/$sample_name/raw_read_QC_reports", mode: 'copy'

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

    stub:
    """
    touch ${sample_name}_fastqc.html
    touch ${sample_name}_fastqc.zip
    """
}

process kraken2 {
    /**
    * @QCcheckpoint only pass to Mykrobe if Kraken's top family classification is Mycobacteriaceae
    */

    tag { sample_name }
    label 'preprocessing'
    label 'normal_cpu'
    label 'high_memory'

    publishDir "${params.output_dir}/$sample_name/output_reads", mode: 'copy', pattern: '*.fq.gz', overwrite: 'true'
    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP", mode: 'copy', pattern: '*_kraken_report.*'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*_err.json'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_reads)
    path(database)

    when:
    enough_reads =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_kraken_report.txt"), path("${sample_name}_kraken_report.json"), emit: kraken2_report
    tuple val(sample_name), path("${sample_name}_cleaned_1.fq.gz"), path("${sample_name}_cleaned_2.fq.gz"), stdout, emit: kraken2_fqs
    path("${sample_name}.err", emit: kraken2_log)

    script:
    kraken2_report = "${sample_name}_kraken_report.txt"
    kraken2_json = "${sample_name}_kraken_report.json"
    kraken2_read_classification = "${sample_name}_read_classifications.txt"
    nonBac_depleted_reads_1 = "${sample_name}_cleaned_1.fq"
    nonBac_depleted_reads_2 = "${sample_name}_cleaned_2.fq"
    error_log = "${sample_name}_err.json"

    """
    kraken2 --threads ${task.cpus} --db . --output ${kraken2_read_classification} --report ${kraken2_report} --paired $fq1 $fq2

    python3 ${baseDir}/bin/parse_kraken_report2.py ${kraken2_report} ${kraken2_json} 0.5 5000

    ${baseDir}/bin/extract_kraken_reads.py -k ${kraken2_read_classification} -r ${kraken2_report} -s $fq1 -s2 $fq2 -o ${nonBac_depleted_reads_1} -o2 ${nonBac_depleted_reads_2} --taxid 2 --include-children --fastq-output >/dev/null

    gzip -f ${nonBac_depleted_reads_1}
    gzip -f ${nonBac_depleted_reads_2}

    rm -rf ${sample_name}_read_classifications.txt

    run_mykrobe=\$(jq '.Mykrobe' ${kraken2_json})

    if [ \$run_mykrobe == '\"true\"' ]; then printf "" >> ${error_log} && printf "${sample_name}"; else echo '{"error":"Kraken's top family hit either wasn't Mycobacteriaceae, or there were < 100k Mycobacteriaceae reads"}' | jq '.' >> ${error_log} && printf "no"; fi
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
  // label 'preprocessing'
  label 'normal_cpu'
  label 'medium_memory'

  publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP", mode: 'copy', pattern: '*.json'

  input:
  tuple val(sample_name), path(fq1), path(fq2), val(run_afanc)
  path(afanc_myco_db)

  when:
  run_afanc =~ /${sample_name}/

  output:
  tuple val(sample_name), path("${sample_name}_afanc_report.json"), stdout, emit: afanc_report
  // tuple val(sample_name), path(fq1), path(fq2), stdout, emit: afanc_fqs

  script:
  afanc_report = "${sample_name}_afanc_report.json"

  """
  afanc screen ${afanc_myco_db} ${fq1} ${fq2} -p 5.0 -n 1000 -o ${sample_name} -t ${task.cpus}
  python3 ${baseDir}/bin/reformat_afanc_json.py ${sample_name}/${sample_name}.json
  printf ${sample_name}
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

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(run_mykrobe)

    when:
    run_mykrobe =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_mykrobe_report.json"), stdout, emit: mykrobe_report
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: mykrobe_fqs

    script:
    mykrobe_report = "${sample_name}_mykrobe_report.json"

    """
    mykrobe predict --sample ${sample_name} --species tb --threads ${task.cpus} --format json --output ${mykrobe_report} -1 $fq1 $fq2
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
    label 'low_memory'

    publishDir "${params.output_dir}/$sample_name/output_reads", mode: 'copy', pattern: '*.fq.gz', overwrite: 'true'

    input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_myco_reads)
    path(index)

    when:
    enough_myco_reads =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_cleaned_1.fq.gz"), path("${sample_name}_cleaned_2.fq.gz"), emit: bowtie2_fqs

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

    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP", mode: 'copy', pattern: '*.json'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*_err.json'

    input:
    tuple val(sample_name), path(fq1), path(fq2), path(mykrobe_json), val(enough_myco_reads), path(kraken_report), path(kraken_json)

    when:
    enough_myco_reads =~ /${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_urllist.txt"), stdout, emit: contam_list
    tuple val(sample_name), path("${sample_name}_species_in_sample_previous.json"), stdout, emit: prev_sample_json
    tuple val(sample_name), path("${sample_name}_species_in_sample.json"), stdout, emit: sample_json
    tuple val(sample_name), path("${sample_name}_nocontam_1.fq.gz"), path("${sample_name}_nocontam_2.fq.gz"), path("${sample_name}_species_in_sample.json"), stdout, emit: nocontam_fqs optional true
    path("${sample_name}_err.json", emit: contam_log)

    script:
    error_log = "${sample_name}_err.json"

    """
    python3 ${baseDir}/bin/identify_tophit_and_contaminants2.py ${mykrobe_json} ${kraken_json} ${baseDir}/resources/assembly_summary_refseq.txt ${params.species} ${params.unmix_myco} ${baseDir}/resources null

    cp ${sample_name}_species_in_sample.json ${sample_name}_species_in_sample_previous.json

    contam_to_remove=\$(jq -r '.summary_questions.are_there_contaminants' ${sample_name}_species_in_sample.json)
    acceptable_species=\$(jq -r '.summary_questions.is_the_top_species_appropriate' ${sample_name}_species_in_sample.json)
    top_hit=\$(jq -r '.top_hit.name' ${sample_name}_species_in_sample.json)

    if [ \$contam_to_remove == 'yes' ]; then printf "NOW_DECONTAMINATE_${sample_name}" && printf "" >> ${error_log}; elif [ \$contam_to_remove == 'no' ] && [ \$acceptable_species == 'yes' ]; then printf "NOW_ALIGN_TO_REF_${sample_name}" && mv $fq1 ${sample_name}_nocontam_1.fq.gz && mv $fq2 ${sample_name}_nocontam_2.fq.gz && printf "" >> ${error_log}; elif [ \$contam_to_remove == 'no' ] && [ \$acceptable_species == 'no' ]; then echo '{"error":"top hit (\$top_hit) is not one of the 10 accepted mycobacteria"}' | jq '.' >> ${error_log}; fi
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

    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*_err.json'

    input:
    tuple val(sample_name), path(contam_list), val(run_decontaminator)

    when:
    run_decontaminator =~ /NOW\_DECONTAMINATE\_${sample_name}/

    output:
    tuple val(sample_name), path("${sample_name}_contaminants.fa"), stdout, emit: contam_fa
    path("${sample_name}.err", emit: downcontam_log)

    script:
    contaminant_fa = "${sample_name}_contaminants.fa"
    error_log = "${sample_name}_err.json"

    """
    wget -i ${contam_list} --spider -nv -a linktestlog.txt 2>&1
    grep -o 'ftp://.*fna.gz' linktestlog.txt > confirmedurllist.txt

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

    if (( \$num_urls_in == \$num_urls_out )); then printf "" >> ${error_log} && printf "${sample_name}"; else echo '{"error":"there were \$num_urls_in contaminant genomes but only \$num_urls_out could be downloaded"}' | jq '.' >> ${error_log} && printf "fail"; fi
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
    tuple val(sample_name), path(fq1), path(fq2), path(contam_fa), val(does_fa_pass)

    when:
    does_fa_pass =~ /${sample_name}/

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

    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP_and_postContamRemoval", mode: 'copy', pattern: '*_kraken_report.*'

    input:
    tuple val(sample_name), path(fq1), path(fq2)
    path(database)

    output:
    tuple val(sample_name), path("${sample_name}_kraken_report.txt"), path("${sample_name}_kraken_report.json"), emit: reKraken_report

    script:
    kraken2_report = "${sample_name}_kraken_report.txt"
    kraken2_json = "${sample_name}_kraken_report.json"
    kraken2_read_classification = "${sample_name}_read_classifications.txt"

    """
    kraken2 --threads ${task.cpus} --db . --output ${kraken2_read_classification} --report ${kraken2_report} --paired $fq1 $fq2

    python3 ${baseDir}/bin/parse_kraken_report2.py ${kraken2_report} ${kraken2_json} 0.5 5000
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
  // label 'preprocessing'
  label 'normal_cpu'
  label 'medium_memory'

  publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP_and_postContamRemoval", mode: 'copy', pattern: '*.json'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(afanc_myco_db)

  output:
  tuple val(sample_name), path("${sample_name}_afanc_report.json"), emit: reAfanc_report
  // tuple val(sample_name), path(fq1), path(fq2), stdout, emit: afanc_fqs

  script:
  afanc_report = "${sample_name}_afanc_report.json"

  """
  afanc screen ${afanc_myco_db} ${fq1} ${fq2} -p 5.0 -n 1000 -o ${sample_name} -t ${task.cpus}
  python3 ${baseDir}/bin/reformat_afanc_json.py ${sample_name}/${sample_name}.json
  printf ${sample_name}
  """

  stub:
  afanc_report = "${sample_name}_afanc_report.json"

  """
  touch ${afanc_report}
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

    input:
    tuple val(sample_name), path(fq1), path(fq2)

    output:
    tuple val(sample_name), path("${sample_name}_mykrobe_report.json"), emit: reMykrobe_report

    script:
    mykrobe_report = "${sample_name}_mykrobe_report.json"

    """
    mykrobe predict --sample ${sample_name} --species tb --threads ${task.cpus} --format json --output ${mykrobe_report} -1 $fq1 $fq2
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

    publishDir "${params.output_dir}/$sample_name/speciation_reports_for_reads_postFastP_and_postContamRemoval", mode: 'copy', pattern: '*.json'
    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*.err'

    input:
    tuple val(sample_name), path(mykrobe_json), path(kraken_report), path(kraken_json), path(prev_species_json), val(decontam)

    output:
    tuple val(sample_name), path("${sample_name}_species_in_sample.json"), stdout, emit: summary_json
    path("${sample_name}.err", emit: summary_log)

    script:
    error_log = "${sample_name}.err"

    """
    python3 ${baseDir}/bin/identify_tophit_and_contaminants2.py ${mykrobe_json} ${kraken_json} ${baseDir}/resources/assembly_summary_refseq.txt ${params.species} ${params.unmix_myco} ${baseDir}/resources ${prev_species_json}

    contam_to_remove=\$(jq -r '.summary_questions.are_there_contaminants' ${sample_name}_species_in_sample.json)
    acceptable_species=\$(jq -r '.summary_questions.is_the_top_species_appropriate' ${sample_name}_species_in_sample.json)
    top_hit=\$(jq -r '.top_hit.name' ${sample_name}_species_in_sample.json)

    if [ \$contam_to_remove == 'yes' ]; then echo '{"error":"sample remains contaminated, even after attempting to resolve this"}' | jq '.' >> ${error_log}; fi

    if [ \$contam_to_remove == 'no' ] && [ \$acceptable_species == 'yes' ]; then printf "NOW_ALIGN_TO_REF_${sample_name}" && printf "" >> ${error_log}; elif [ \$contam_to_remove == 'no' ] && [ \$acceptable_species == 'no' ]; then echo '{"error":"top hit (\$top_hit) is not one of the 10 accepted mycobacteria"}' | jq '.' >> ${error_log}; fi
    """

    stub:
    error_log = "${sample_name}.err"

    """
    touch ${sample_name}_species_in_sample.json
    touch ${error_log}
    printf ${params.summary_doWeAlign}
    """
}
