// enable dsl2
nextflow.enable.dsl = 2

include {identifyBacterialContaminants} from '../modules/decontaminationModules.nf' params(params)
include {downloadContamGenomes} from '../modules/decontaminationModules.nf' params(params)
include {mapToContamFa} from '../modules/decontaminationModules.nf' params(params)
include {reKraken} from '../modules/decontaminationModules.nf' params(params)
include {reAfanc} from '../modules/decontaminationModules.nf' params(params)
include {reMykrobe} from '../modules/decontaminationModules.nf' params(params)
include {summarise} from '../modules/decontaminationModules.nf' params(params)
include {count_pass} from '../modules/decontaminationModules.nf' params(params)

workflow decontaminate {
    take:
    fastqs // tuple containing sample name, fq1 and fq2
    krakenDB //kraken datatase, won't change
    afanc_myco_db //afanc database, won't change
    speciation_report //afanc report, recursion takes from reafanc
    kraken2_json //kraken report, recursion takes from rekraken
    resource_dir //directory with references etc. won't change 
    refseq_path //path to refseq file won't change 

    //tracking variables
    pass_number //integeter of number of iterations
    do_we_decontaminate //bool, have we finished?

    main:
      
      //sample name, fastq1, fastq2, software versions, afanc_report, do_we_attempt, kraken text, kraken json
      contaminant_input = fastqs.join(speciation_report, by: 0).join(kraken2_json, by: 0)
      
      //just sample name
      sample_name = fastqs.map{it[0]}

      //Identify contaminants from kraken report, download the contaminants and align against bwa
      identifyBacterialContaminants(contaminant_input, resource_dir, refseq_path, pass_number)

      
      //overwrite do_we_break and summary_json incase we break early (e.g. there are no contams)
      summary_json = identifyBacterialContaminants.out.nocontam_fqs.map{it[4]}
      do_we_break  = identifyBacterialContaminants.out.nocontam_fqs.map{it[5]}

      downloadContamGenomes(identifyBacterialContaminants.out.contam_list, pass_number)
      mapToContamFa(fastqs.join(downloadContamGenomes.out.contam_fa, by: 0), pass_number)
        
      //override fastqs parameter with our new reads
      fastqs = mapToContamFa.out.reClassification_fqs

      //re-run our speciation tools 
      reKraken(fastqs, krakenDB.toList(), pass_number)
      reMykrobe(fastqs, pass_number)
      reAfanc(fastqs, afanc_myco_db, pass_number)
      
      //set speciation report
      speciation_report = reAfanc.out.reAfanc_report
      kraken_report     = reKraken.out.reKraken_report
      previous_json     = identifyBacterialContaminants.out.prev_sample_json

      //parse reports and make summary json
      summarise(speciation_report.join(kraken_report, by: 0).join(previous_json, by: 0), resource_dir, refseq_path, pass_number)
      do_we_break = summarise.out.do_we_break


    emit:
      decontamination_json  = summarise.out.summary_json
      do_we_proceed         = do_we_break
      decontaminated_fastqs = mapToContamFa.out.reClassification_fqs
      speciation_report     = speciation_report
      kraken_report         = kraken_report
      previous_json         = previous_json
      uncontaminatd_fastqs  = identifyBacterialContaminants.out.nocontam_fqs
}
