// enable dsl2
nextflow.enable.dsl = 2

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
      }

      if ( params.filetype == "fastq" ) {

          checkFqValidity(input_files)

          countReads(checkFqValidity.out.checkValidity_fqs)
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

      decontaminate(bowtie2.out.bowtie2_fqs, krakenDB, afanc_myco_db, bowtie2.out.bowtie2_fqs, speciation_report,
                    kraken2.out.kraken2_json, resource_dir, refseq_path, bowtie2.out.bowtie2_fqs, kraken2.out.kraken2_json, )


    emit:
    
      decontam_seqs = decontaminate.out.reClassification_fqs
      decontam_json = decontaminate.out.summary_json
      nocontam_seqs_json = decontaminate.out.nocontam_fqs
}
