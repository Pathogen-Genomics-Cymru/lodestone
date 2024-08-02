// enable dsl2
nextflow.enable.dsl = 2
nextflow.preview.recursion=true

// import modules
include {checkFqValidity} from '../modules/preprocessingModules.nf' params(params)
include {countReads} from '../modules/preprocessingModules.nf' params(params)
include {fastp} from '../modules/preprocessingModules.nf' params(params)
include {fastQC} from '../modules/preprocessingModules.nf' params(params)
include {kraken2} from '../modules/preprocessingModules.nf' params(params)
include {afanc} from '../modules/preprocessingModules.nf' params(params)
include {mykrobe} from '../modules/preprocessingModules.nf' params(params)
include {bowtie2} from '../modules/preprocessingModules.nf' params(params)
include {checkBamValidity} from '../modules/preprocessingModules.nf' params(params)
include {bam2fastq} from '../modules/preprocessingModules.nf' params(params)

//import subworkflow
include {decontaminate} from '../subworkflows/decontamination.nf'

// define workflow component
workflow preprocessing {

    take:
      input_files
      krakenDB
      bowtie_dir
      afanc_myco_db
      resource_dir
      refseq_path


    main:

      if ( params.filetype == "bam" ) {

          checkBamValidity(input_files)

          bam2fastq(checkBamValidity.out.checkValidity_bam)

          countReads(bam2fastq.out.bam2fastq_fqs)

          versions = input_files.map{it[2]}
      }

      if ( params.filetype == "fastq" ) {

          checkFqValidity(input_files)

          countReads(checkFqValidity.out.checkValidity_fqs)

          versions = input_files.map{it[3]}
      }

      fastp(countReads.out.countReads_fqs)

      fastQC(fastp.out.fastp_fqs)

      kraken2(fastp.out.fastp_fqs, krakenDB.toList())

      mykrobe(kraken2.out.kraken2_fqs)
      
      afanc(kraken2.out.kraken2_fqs.join(kraken2.out.kraken2_json, by: 0), afanc_myco_db, resource_dir, refseq_path)
      
      // set speciation report
      speciation_report = afanc.out.afanc_json

      bowtie2(kraken2.out.kraken2_fqs, bowtie_dir.toList())

      //subworkflow to remove erraneous species
      pass_number = 1
      decontamination_finshed = "FALSE"
      
      decontaminate.recurse(bowtie2.out.bowtie2_fqs.collect(), bowtie2.out.software_json.collect(), krakenDB.collect(), afanc_myco_db, 
                            speciation_report.collect(), kraken2.out.kraken2_json.collect(), resource_dir, refseq_path, 
                            versions.collect(), pass_number, decontamination_finshed).until{it[10] != "FALSE"}

      /*to pass to clockwork we need the following:
        - sample name
        - fq1 (cleaned or not)
        - fq2
        - versions json
        - species in sample (e.g. summarise report)
        - do we align bool
     */

      sample_and_fastq = decontaminate.out.fastqs
      versions = versions
      summary_json = decontaminate.out.summary_json.map{it[1]} //we've kept the samee summary json for the recursion, so now we can just grab the json
      do_we_align = decontaminate.out.decontamination_finshed
      output_to_pass = sample_and_fastq.combine(versions).combine(summary_json).combine(do_we_align)

    emit:
    fastqs_and_reports = output_to_pass
}
