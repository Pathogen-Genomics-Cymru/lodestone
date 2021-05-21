// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {checkFqValidity} from '../modules/preprocessingModules.nf' params(params)
include {countReads} from '../modules/preprocessingModules.nf' params(params)
include {fastp} from '../modules/preprocessingModules.nf' params(params)
include {fastQC} from '../modules/preprocessingModules.nf' params(params)
include {kraken2} from '../modules/preprocessingModules.nf' params(params)
include {mykrobe} from '../modules/preprocessingModules.nf' params(params)
include {bowtie2} from '../modules/preprocessingModules.nf' params(params)
include {identifyBacterialContaminants} from '../modules/preprocessingModules.nf' params(params)
include {downloadContamGenomes} from '../modules/preprocessingModules.nf' params(params)
include {mapToContamFa} from '../modules/preprocessingModules.nf' params(params)
include {reKraken} from '../modules/preprocessingModules.nf' params(params)
include {reMykrobe} from '../modules/preprocessingModules.nf' params(params)
include {summarise} from '../modules/preprocessingModules.nf' params(params)
include {checkBamValidity} from '../modules/preprocessingModules.nf' params(params)
include {bam2fastq} from '../modules/preprocessingModules.nf' params(params)

// define workflow component
workflow preprocessing {

    take:
      input_files
      krakenDB
      bowtie_dir

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

      bowtie2(kraken2.out.kraken2_fqs, bowtie_dir)

      identifyBacterialContaminants(mykrobe.out.mykrobe_report, kraken2.out.kraken2_report)

      downloadContamGenomes(identifyBacterialContaminants.out.contam_list)

      mapToContamFa(bowtie2.out.bowtie2_fqs, downloadContamGenomes.out.contam_fa)

      reKraken(mapToContamFa.out.reClassification_fqs, krakenDB.toList())

      reMykrobe(mapToContamFa.out.reClassification_fqs)

      summarise(reMykrobe.out.reMykrobe_report, reKraken.out.reKraken_report, identifyBacterialContaminants.out.prev_sample_json)
}
