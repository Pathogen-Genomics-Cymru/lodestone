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

// define workflow component
workflow vcfpredict {

    take:
      sample_and_fastqs
      clockwork_bam
      clockwork_bcftools_tuple
      minos_vcf_tuple
      reference_fasta
      

    main:

      if ( params.vcfmix == "yes" ) {

          vcfmix(clockwork_bcftools_tuple)

      }

      //get just the vcf
      sample_name = minos_vcf_tuple.map{it[0]}
      minos_vcf = minos_vcf_tuple.map{it[1]}
      do_we_resistance_profile = minos_vcf_tuple.map{it[2]}
      report_json  = minos_vcf_tuple.map{it[3]}
      bam = clockwork_bam.map{it[2]}
      fastq_and_report = sample_and_fastqs.combine(report_json).combine(do_we_resistance_profile)

      if ( params.resistance_profiler == "tb-profiler"){

        //if we are local and want to match our references, run this
        if (params.update_tbprofiler == "yes"){
        tbprofiler_update_db(reference_fasta)
        }
        
        //add allelic depth back in: was calculated in mpileup but lost in minos
        add_allelic_depth(sample_name, minos_vcf, bam, reference_fasta, do_we_resistance_profile)
        //run tb-profiler
        tbprofiler(sample_name, add_allelic_depth.out, report_json, do_we_resistance_profile)
        profiling_json = tbprofiler.out.tbprofiler_json
        if(params.collate == "yes"){
          collated_jsons = tbtamr.out.collate_json.collect()
          tbprofiler_collate(collated_jsons)
        }
      } else if (params.resistance_profiler == "tbtamr"){
        tbtamr(fastq_and_report)
        profiling_json = tbtamr.out.tbtamr_json
        if(params.collate == "yes"){
          collated_jsons = tbtamr.out.collate_json.collect()
          tbtamr_collate(collated_jsons)
        }
      }
      
      if (params.vcfmix == "yes" && params.resistance_profiler != "none"){
          finalJson(vcfmix.out.vcfmix_json.join(profiling_json, by: 0))
      }
}
