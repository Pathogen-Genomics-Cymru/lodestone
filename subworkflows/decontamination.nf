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
    bowtie_fastqs
    krakenDB
    afanc_myco_db
    bowtie_fqs
    speciation_report
    kraken2_json
    resource_dir
    refseq_path

    //downstream inputs (e.g. what will get outputted later)
    reClassification_fqs
    summary_json
    nocontam_fqs

    //tracking variables
    pass_number
    decontamination_finshed

    main:
      identifyBacterialContaminants(bowtie_fqs.join(speciation_report, by: 0).join(kraken2_json, by: 0), resource_dir, refseq_path)
      downloadContamGenomes(identifyBacterialContaminants.out.contam_list)
      mapToContamFa(bowtie_fastqs.join(downloadContamGenomes.out.contam_fa, by: 0))
      
      //re-apply kraken, mykrobe and afanc based on our (hopefully now clean) reads
      reKraken(mapToContamFa.out.reClassification_fqs, krakenDB.toList())
      reMykrobe(mapToContamFa.out.reClassification_fqs)
      reAfanc(mapToContamFa.out.reClassification_fqs, afanc_myco_db)

      // set speciation report
      speciation_report = reAfanc.out.reAfanc_report
      
      // parse reports and make summary json
      summarise(speciation_report.join(reKraken.out.reKraken_report, by: 0).join(identifyBacterialContaminants.out.prev_sample_json, by: 0), resource_dir, refseq_path)

      new_pass_number = pass_number + 1

    emit:
    afanc_myco_db = afanc_myco_db
    bowtie_fqs = mapToContamFa.out.reClassification_fqs
    speciation_report = reAfanc.out.afanc_json
    kraken2_json = reKraken.out.kraken2_json
    resource_dir = resource_dir
    refseq_path = refseq_path

    //real output
    reClassification_fqs = mapToContamFa.out.reClassification_fqs
    summary_json = summarise.out.summary_json
    nocontam_fqs = identifyBacterialContaminants.out.nocontam_fqs

    //tracking variables
    pass_number = new_pass_number
    decontamination_finshed = do_we_break

}