#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

//nf-validation
include { validateParameters; paramsHelp; paramsSummaryLog} from 'plugin/nf-schema'

// import subworkflows
include {preprocessing} from './workflows/preprocessing.nf'
include {clockwork} from './workflows/clockwork.nf'
include {vcfpredict} from './workflows/vcfpredict.nf'
include {getversion} from './workflows/getversion.nf'

def helpMessage() {
log.info paramsHelp('Nextflow run main.nf -profile singularity --filetype fastq --input_dir fq_dir --pattern "*_R{1,2}.fastq.gz" --unmix_myco yes --output_dir ./')
}

if (params.help) {
    helpMessage()
    exit(0)
}

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

// main workflow
workflow {

    // add a trailing slash if it was not originally provided to --input_dir
    inputdir_amended = "${params.input_dir}".replaceFirst(/$/, "/")

    indir = inputdir_amended
    numfiles = 0

    if ( params.filetype == "bam" ) {
        reads = indir + "*.bam"
        numfiles = file(reads) // count the number of files

        Channel.fromPath(reads)
               .set{ input_files }
    }

    if ( params.filetype == "fastq" ) {
        pattern = params.pattern
	reads = indir + pattern
	numfiles = file(reads) // count the number of files

        Channel.fromFilePairs(reads, flat: true, checkIfExists: true, size: -1)
	       .ifEmpty { error "cannot find any reads matching ${pattern} in ${indir}" }
	       .set{ input_files }
    }

    // create channels for kraken2 database and bowtie2 index
    krakenDB = Channel.fromPath( "${params.kraken.kraken_db}/*.k2d" )
    bowtie_dir = Channel.fromPath( "${params.bowtie.bowtie_index}/*.bt2" )

    // main workflow
    main:

      // GETVERSION SUB-WORKFLOW

      getversion()

      // PREPROCESSING SUB-WORKFLOW

      input_files_vjson = input_files.combine(getversion.out.getversion_json)

      preprocessing(input_files_vjson, krakenDB, bowtie_dir, params.afanc.afanc_myco_db, 
                    params.resources.resource_dir, params.resources.refseq)

      // CLOCKWORK SUB-WORKFLOW
      preprocessing_output = preprocessing.out.fastqs_and_reports
      clockwork(preprocessing_output)

      // VCFPREDICT SUB-WORKFLOW
      profiler_input_vcf = clockwork.out.profiler_input_vcf 
      profiler_input_fq = clockwork.out.profiler_input_fq     

      vcfpredict(profiler_input_fq, profiler_input_vcf)
    
}

workflow.onComplete {
    if ( workflow.success ) {
        log.info """
        ===========================================
        ${ANSI_GREEN}Workflow completed successfully
        """
        .stripIndent()
    }
    else {
        log.info """
        ===========================================
        ${ANSI_RED}Finished with errors${ANSI_RESET}
        """
        .stripIndent()
    }
}
