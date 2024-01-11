// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {vcfmix} from '../modules/vcfpredictModules.nf' params(params)
include {tbprofiler} from '../modules/vcfpredictModules.nf' params(params)
include {tbprofiler_update_db} from '../modules/vcfpredictModules.nf' params(params)
include {add_allelic_depth} from '../modules/vcfpredictModules.nf' params(params) 

// define workflow component
workflow vcfpredict {

    take:
      clockwork_bcftools_tuple
      minos_vcf_tuple
      reference_fasta
      

    main:

      if ( params.vcfmix == "yes" ) {

          vcfmix(clockwork_bcftools_tuple)

      }

      if ( params.resistance_profiler == "tb-profiler"){
        //get just the vcf
        sample_name = minos_vcf_tuple.map{it[0]}
        minos_vcf = minos_vcf_tuple.map{it[1]}
        do_we_resistance_profile = minos_vcf_tuple.map{it[2]}
        report_json  = minos_vcf_tuple.map{it[3]}

        if (params.update_tbprofiler == "yes"){
        tbprofiler_update_db(reference_fasta)
        }
        
        //add allelic depth back in: was calculated in mpileup but lost in minos
        add_allelic_depth(sample_name, minos_vcf, reference_fasta, do_we_resistance_profile)
        tbprofiler(sample_name, add_allelic_depth,out, report_json, do_we_resistance_profile)
      }
      
      if (params.vcfmix == "yes" && params.resistance_profiler != "none"){
          finalJson(vcfmix.out.vcfmix_json.join(gnomonicus.out.tbprofiler_json, by: 0))
      }
}
