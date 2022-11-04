// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {vcfmix} from '../modules/vcfpredictModules.nf' params(params)
include {gnomon} from '../modules/vcfpredictModules.nf' params(params)

// define workflow component
workflow vcfpredict {

    take:

      clockwork_bcftools
      clockwork_minos
      resources_dir

    main:

      if ( params.vcfmix == "yes" ) {

          vcfmix(clockwork_bcftools)

      }

      if ( params.gnomon == "yes" ) {

          gnomon(clockwork_minos, resources_dir)

      }

}
