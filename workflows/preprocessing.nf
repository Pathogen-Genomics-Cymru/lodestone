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
include {decontaminate as decontaminate_second_pass} from '../subworkflows/decontamination.nf'

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
      decontamination_finished = "FALSE"
      software_json = bowtie2.out.software_json

      first_pass_decontaminate = decontaminate(bowtie2.out.bowtie2_fqs, krakenDB, afanc_myco_db, 
                    speciation_report, kraken2.out.kraken2_json, resource_dir, 
                    refseq_path, 1, decontamination_finished)


      /*to pass to clockwork we need the following:
        - sample name
        - fq1 (cleaned or not)
        - fq2
        - versions json
        - species in sample (e.g. summarise report)
        - do we align bool
        write them first
       */

      fastqs          = first_pass_decontaminate.decontaminated_fastqs
      summary_json    = first_pass_decontaminate.decontamination_json
      do_we_repass    = first_pass_decontaminate.do_we_proceed
      nomix_seqs_json = first_pass_decontaminate.uncontaminatd_fastqs

      //grab new reports for our second iteration
      //replace original sample_id val that was used to start identifyBacterialContaminants with our last pass
      reafanc_report  = first_pass_decontaminate.speciation_report.combine(do_we_repass) 
      rekraken_report = first_pass_decontaminate.kraken_report

      first_pass_output = fastqs.join(summary_json, by: 0).mix(nomix_seqs_json)

      second_pass_decontaminate = decontaminate_second_pass(fastqs, krakenDB, afanc_myco_db, 
                                                            reafanc_report, rekraken_report, resource_dir, 
                                                            refseq_path, 2, decontamination_finished)

      fastqs_2nd          = second_pass_decontaminate.decontaminated_fastqs
      summary_json_2nd    = second_pass_decontaminate.decontamination_json
      do_we_align         = second_pass_decontaminate.do_we_proceed        
      nomix_seqs_json_2nd = second_pass_decontaminate.uncontaminatd_fastqs

      //combine to suit how clockwork wf was written
      second_pass_output = fastqs_2nd.join(summary_json_2nd, by: 0).mix(nomix_seqs_json_2nd)

      both_passes = first_pass_output.concat(second_pass_output)

    emit:
    fastqs_and_reports = both_passes
}
