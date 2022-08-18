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

    main:

      vcfmix(clockwork_bcftools)

      gnomon(clockwork_minos)

}
