// enable dsl2
nextflow.enable.dsl = 2

include {identifyBacterialContaminants} from '../modules/preprocessingModules.nf' params(params)
include {downloadContamGenomes} from '../modules/preprocessingModules.nf' params(params)
include {mapToContamFa} from '../modules/preprocessingModules.nf' params(params)
include {reKraken} from '../modules/preprocessingModules.nf' params(params)
include {reAfanc} from '../modules/preprocessingModules.nf' params(params)
include {reMykrobe} from '../modules/preprocessingModules.nf' params(params)
include {summarise} from '../modules/preprocessingModules.nf' params(params)

//because this script is getting ran recursively, input and output need to be the same

workflow decontaminate {
    take:
    //upstream inputs (e.g. what we need to actually start process)
    fastqs // tuple containing sample name, fq1 and fq2
    summary_json //sample json from bowtie, will take from summarise
    krakenDB //kraken datatase, won't change
    afanc_myco_db //afanc database, won't change
    speciation_report //afanc report, recursion takes from reafanc
    kraken2_json //kraken report, recursion takes from rekraken
    resource_dir //directory with references etc. won't change 
    refseq_path //path to refseq file won't change 
    version_json //version json file from first function.. won't change

    //tracking variables
    pass_number //integeter of number of iterations
    decontamination_finshed //bool, have we finished?

    main:
      fastq_and_version = fastqs.combine(version_json)
      contaminant_input = fastq_and_version.join(speciation_report, by: 0).join(kraken2_json, by: 0) //tuple val(sample_name), path(fq1), path(fq2), path(software_json), 
                                                                                                    //path(afanc_json), val(enough_myco_reads), path(kraken_report), path(kraken_json)
      
      identifyBacterialContaminants(contaminant_input, resource_dir, refseq_path, pass_number)
      downloadContamGenomes(identifyBacterialContaminants.out.contam_list, pass_number)
      mapToContamFa(fastq_and_version.join(downloadContamGenomes.out.contam_fa, by: 0), pass_number)
      
      //re-apply kraken, mykrobe and afanc based on our (hopefully now clean) reads
      
      //override fastqs parameter with our new reads
      fastqs = mapToContamFa.out.reClassification_fqs
      fastq_and_version = fastqs.combine(version_json)

      reKraken(fastq_and_version, krakenDB, pass_number)
      reMykrobe(fastq_and_version, pass_number)
      reAfanc(fastq_and_version, afanc_myco_db, pass_number)

      // set speciation report
      speciation_report = reAfanc.out.reAfanc_report
      
      // parse reports and make summary json
      summarise(speciation_report.join(reKraken.out.reKraken_report, by: 0).join(identifyBacterialContaminants.out.prev_sample_json, by: 0), resource_dir, refseq_path, pass_number)

      //a rather convoluted way of tracking the decontamination passes because just adding doesn't want to work.. Can perhaps use $task.index
      new_pass_number = pass_number.combine(channel.of(1)).sum()

    emit:
    fastqs
    summary_json = summarise.out.summary_json
    krakenDB
    afanc_myco_db
    speciation_report = reAfanc.out.reAfanc_report
    kraken2_json = reKraken.out.reKraken_report
    resource_dir
    refseq_path
    version_json 

    pass_number = new_pass_number
    decontamination_finshed = summarise.out.do_we_break
}