// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {vcfmix} from '../modules/vcfpredictModules.nf' params(params)
include {tbprofiler} from '../modules/vcfpredictModules.nf' params(params)
include {tbprofiler_update_db} from '../modules/vcfpredictModules.nf' params(params)
include {add_allelic_depth} from '../modules/vcfpredictModules.nf' params(params) 
include {finalJson} from '../modules/vcfpredictModules.nf' params(params) 
include {tbtamr} from '../modules/vcfpredictModules.nf' params(params)
include {tbtamr_collate} from '../modules/vcfpredictModules.nf' params(params)
include {tbprofiler_collate} from '../modules/vcfpredictModules.nf' params(params)
include {ntmprofiler} from '../modules/vcfpredictModules.nf' params(params)
include {ntmprofiler_collate} from '../modules/vcfpredictModules.nf' params(params)
// define workflow component
workflow vcfpredict {

    take:
    profiler_input_fq
    profiler_input_vcf

    main:
      //ntm-profiling: e.g. everything down being passed into tbtamr/tb-profiler
      //at the moment it is only ran on fastqs; need to find a sensible way
      //of linking up the references
      ntmprofiler(profiler_input_fq)

      ntm_profiling_out = ntmprofiler.out.vcfmix_in
      
      if(params.collate == "yes"){
        collated_ntm_jsons = ntmprofiler.out.collate_json.collect()
        ntmprofiler_collate(collated_ntm_jsons)
      }

      if ( params.resistance_profiler == "tb-profiler"){

        //if we are local and want to match our references, run this
        if (params.update_tbprofiler == "yes"){
        tbprofiler_update_db(reference_fasta)
        }
        
        //add allelic depth back in: was calculated in mpileup but lost in minos
        add_allelic_depth(profiler_input_vcf)
        //run tb-profiler
        tbprofiler(add_allelic_depth.out)

        tb_profiling_out = tbprofiler.out.vcfmix_in

        if(params.collate == "yes"){
          collated_jsons = tbprofiler.out.collate_json.collect()
          tbprofiler_collate(collated_jsons)
        }
      } else if (params.resistance_profiler == "tbtamr"){
        tbtamr(profiler_input_fq)

        tb_profiling_out = tbtamr.out.vcfmix_in
        
        if(params.collate == "yes"){
          collated_jsons = tbtamr.out.collate_json.collect()
          tbtamr_collate(collated_jsons)
        }
      }
      
      profiling_jsons = ntm_profiling_out.mix(tb_profiling_out)
      profiling_jsons.view()
      vcfmix(profiling_jsons)
}
