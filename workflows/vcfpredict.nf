// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {vcfmix} from '../modules/vcfpredictModules.nf' params(params)
include {gnomonicus} from '../modules/vcfpredictModules.nf' params(params)

// define workflow component
workflow vcfpredict {

    take:

      clockwork_bcftools
      clockwork_minos

    main:

      if ( params.vcfmix == "yes" ) {

          vcfmix(clockwork_bcftools)

      }

      if ( params.gnomonicus == "yes" ) {

          gnomonicus(clockwork_minos)

      }

}
