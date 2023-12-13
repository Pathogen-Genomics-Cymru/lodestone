// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {vcfmix} from '../modules/vcfpredictModules.nf' params(params)
include {tbprofiler} from '../modules/vcfpredictModules.nf' params(params)
include {tbprofiler_update_db} from '../modules/vcfpredictModules.nf' params(params)

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
        minos_vcf = minos_vcf_tuple.map{it[1]}
        sample_name = minos_vcf_tuple.map{it[0]}

        tbprofiler_update_db(reference_fasta)
        tbprofiler(sample_name, minos_vcf)
      }
}
