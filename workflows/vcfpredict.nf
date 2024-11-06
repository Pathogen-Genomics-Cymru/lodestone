// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {vcfmix} from '../modules/vcfpredictModules.nf'
include {tbprofiler} from '../modules/vcfpredictModules.nf'
include {tbprofiler_update_db} from '../modules/vcfpredictModules.nf'
include {add_allelic_depth} from '../modules/vcfpredictModules.nf' 
include {finalJson} from '../modules/vcfpredictModules.nf' 
include {tbtamr} from '../modules/vcfpredictModules.nf'
include {tbtamr_collate} from '../modules/vcfpredictModules.nf'
include {tbprofiler_collate} from '../modules/vcfpredictModules.nf'
include {ntmprofiler} from '../modules/vcfpredictModules.nf'
include {ntmprofiler_collate} from '../modules/vcfpredictModules.nf'
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
      
      if(params.resistance.collate == true){
        collated_ntm_jsons = ntmprofiler.out.collate_json.collect()
        ntmprofiler_collate(collated_ntm_jsons)
      }

      if ( params.resistance.resistance_profiler == "tb-profiler"){
        
        //run tb-profiler
        tbprofiler(profiler_input_vcf)

        tb_profiling_out = tbprofiler.out.vcfmix_in

        if(params.resistance.collate == true){
          collated_jsons = tbprofiler.out.collate_json.collect()
          tbprofiler_collate(collated_jsons)
        }
      } else if (params.resistance.resistance_profiler == "tbtamr"){
        tbtamr(profiler_input_fq)

        tb_profiling_out = tbtamr.out.vcfmix_in
        
        if(params.resistance.collate == true){
          collated_jsons = tbtamr.out.collate_json.collect()
          tbtamr_collate(collated_jsons)
        }
      }
      
      profiling_jsons = ntm_profiling_out.mix(tb_profiling_out)
      vcfmix(profiling_jsons)
}
