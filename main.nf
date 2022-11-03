#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// import subworkflows
include {preprocessing} from './workflows/preprocessing.nf'
include {clockwork} from './workflows/clockwork.nf'
include {vcfpredict} from './workflows/vcfpredict.nf'
include {getversion} from './workflows/getversion.nf'

/*
 ANSI escape codes to allow colour-coded output messages
 This code is from https://github.com/angelovangel
 */
 
ANSI_GREEN = "\033[1;32m"
ANSI_RED   = "\033[1;31m"
ANSI_RESET = "\033[0m"

if (params.help) {
    helpMessage()
    exit(0)
}

def helpMessage() {
log.info """
========================================================================
M Y C O B A C T E R I A L  P I P E L I N E
  
Cleans and QCs reads with fastp and FastQC, classifies with Kraken2 & Mykrobe, removes non-bacterial content, and - by alignment to any minority genomes - disambiguates mixtures of bacterial reads.

Takes as input one directory containing pairs of fastq(.gz) or bam files.
Produces as output one directory per sample, containing the relevant reports & a pair of cleaned fastqs.
	
Mandatory and conditional parameters:
------------------------------------------------------------------------
--input_dir           Directory containing fastq OR bam files. Workflow will process one or the other, so don't mix
--filetype	      File type in input_dir. One of either "fastq" or "bam". fastq files can be gzipped and do not 
                      have to literally take the form "*.fastq"; see --pattern
--pattern             Regex to match files in input_dir, e.g. "*_R{1,2}.fq.gz". Only mandatory if --filetype is "fastq"
--output_dir          Output directory, in which will be created subdirectories matching base name of fastq/bam files
--unmix_myco	      Do you want to disambiguate mixed-mycobacterial samples by read alignment? One of "yes" or "no"
	              If "yes" workflow will remove reads mapping to any minority mycobacterial genomes but in doing so 
                      WILL ALMOST CERTAINLY ALSO reduce coverage of the principal species
	              If "no" then mixed-mycobacterial samples will be left alone. Mixtures of mycobacteria + non-mycobacteria 
                      will still be disambiguated
--kraken_db           Directory containing Kraken2 database files (obtain from https://benlangmead.github.io/aws-indexes/k2)
--bowtie2_index       Directory containing Bowtie2 index (obtain from ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19_1kgmaj_bt2.zip
                      This is the Langmead lab pre-built major-allele-SNP reference; see https://github.com/BenLangmead/bowtie-majref)
--bowtie_index_name   Name of the bowtie index, e.g. hg19_1kgmaj
--vcfmix	      Run VFCMIX "yes" or "no". Should be set to "no" for synthetic samples
--gnomon              Run gnomon "yes" or "no"
--amr_cat             Path to the AMR catalogue (https://github.com/oxfordmmm/tuberculosis_amr_catalogues is at /tuberculosis_amr_catalogues
                      in the vcfpredict container)
--report_dir          Output directory for execution reports

Optional parameters:
------------------------------------------------------------------------
--species          principal species in each sample, assuming genus Mycobacterium
                   if parameter used, takes 1 of 10 values: abscessus, africanum, avium, bovis, chelonae, chimaera,
                   fortuitum, intracellulare, kansasii, tuberculosis
                   default: null
                   using this parameter will apply an additional sanity test to your sample
				   
	           if you DO NOT use this parameter (default option), pipeline will determine principal species from 
                   the reads and consider any other species a contaminant
                   
	           if you DO use this parameter, pipeline will expect this to be the principal species. It will fail 
		   the sample if reads from this species are not actually the majority

                   
Profiles:
------------------------------------------------------------------------
singularity        to run with singularity
docker             to run with docker
agc                to run with Amazon Genomics CLI

			   
Examples:
------------------------------------------------------------------------
nextflow run main.nf -profile singularity --filetype fastq --input_dir fq_dir --pattern "*_R{1,2}.fastq.gz" --unmix_myco yes --output_dir .
nextflow run main.nf -profile docker --filetype bam --input_dir bam_dir --unmix_myco yes --output_dir .
========================================================================
"""
.stripIndent()
}


// confirm that mandatory parameters have been set and that the conditional parameter, --pattern, has been used appropriately
if ( params.input_dir == "" ) {
    exit 1, "error: --input_dir is mandatory (run with --help to see parameters)"
}
if ( params.filetype == "" ) {
    exit 1, "error: --filetype is mandatory (run with --help to see parameters)"
}
if ( ( params.filetype == "fastq" ) && ( params.pattern == "" ) ) {
    exit 1, "error: --pattern is mandatory if you are providing fastq input; describes files in --input_dir (e.g. \"*_R{1,2}.fastq.gz\") (run with --help to see parameters)"
}
if ( ( params.filetype == "bam" ) && ( params.pattern != "" ) ) {
    exit 1, "error: --pattern should only be set if you are providing fastq input (run with --help to see parameters)"
}
if ( params.output_dir == "" ) {
    exit 1, "error: --output_dir is mandatory (run with --help to see parameters)"
}
if ( ( params.filetype != "fastq" ) && ( params.filetype != "bam" ) ) {
    exit 1, "error: --filetype is mandatory and must be either \"fastq\" or \"bam\""
}
if ( ( params.unmix_myco != "yes" ) && ( params.unmix_myco != "no" ) ) {
    exit 1, "error: --unmix_myco is mandatory and must be either \"yes\" or \"no\""
}
if ( ( params.species != "null" ) && ( params.species != "abscessus" ) && ( params.species != "africanum" ) && ( params.species != "avium" ) && ( params.species != "bovis" ) && ( params.species != "chelonae" ) && ( params.species != "chimaera" ) && ( params.species != "fortuitum" ) && ( params.species != "intracellulare" ) && ( params.species != "kansasii" ) && ( params.species != "tuberculosis" ) ) {
    exit 1, "error: --species is optional, but if used should be one of either abscessus, africanum, avium, bovis, chelonae, chimaera, fortuitum, intracellulare, kansasii, tuberculosis"
}
if ( params.report_dir == "" ) {
    exit 1, "error: --report_dir is mandatory (run with --help to see parameters)"
}

log.info """
========================================================================
M Y C O B A C T E R I A L  P I P E L I N E

Parameters used:
------------------------------------------------------------------------
--input_dir		${params.input_dir}
--filetype		${params.filetype}
--pattern		${params.pattern}
--output_dir	        ${params.output_dir}
--unmix_myco	        ${params.unmix_myco}
--kraken_db		${params.kraken_db}
--bowtie2_index         ${params.bowtie2_index}
--bowtie_index_name     ${params.bowtie_index_name}
--species		${params.species}
--vcfmix		${params.vcfmix}
--gnomon		${params.gnomon}
--amr_cat		${params.amr_cat}
--report_dir		${params.report_dir}

Runtime data:
------------------------------------------------------------------------
Running with profile  ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
Running as user       ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
Launch directory      ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
Workflow session ID   ${ANSI_GREEN}${workflow.sessionId}${ANSI_RESET}
"""
.stripIndent()

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

    // create channel for singularity directory
    sing_dir_ch = Channel.fromPath( "${params.sing_dir}" )

    // create channels for kraken2 database and bowtie2 index
    krakenDB = Channel.fromPath( "${params.kraken_db}/*.k2d" )
    bowtie_dir = Channel.fromPath( "${params.bowtie2_index}/*.bt2" )

    // create channel for resources
    resources_dir_ch = Channel.fromPath( "${params.resources_dir}" )

    // main workflow
    main:

      // GETVERSION SUB-WORKFLOW
      getversion(sing_dir_ch)

      // PREPROCESSING SUB-WORKFLOW
      preprocessing(input_files, krakenDB, bowtie_dir, resources_dir_ch)


      // CLOCKWORK SUB-WORKFLOW
      if ( params.unmix_myco == "yes" ) {

          clockwork_seqs = preprocessing.out.decontam_seqs
          clockwork_json = preprocessing.out.decontam_json

          nomix_seqs_json = preprocessing.out.nocontam_seqs_json

          clockwork(clockwork_seqs.join(clockwork_json, by: 0).mix(nomix_seqs_json), resources_dir_ch)

      }

      if ( params.unmix_myco == "no" ) {

          clockwork_seqs = preprocessing.out.contam_seqs
          clockwork_json = preprocessing.out.contam_json

          clockwork(clockwork_seqs.join(clockwork_json, by: 0))
      }

      mpileup_vcf = clockwork.out.mpileup_vcf
      minos_vcf = clockwork.out.minos_vcf

      // VCFPREDICT SUB-WORKFLOW
      vcfpredict(mpileup_vcf, minos_vcf, resources_dir_ch)

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
