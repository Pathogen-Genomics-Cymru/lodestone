// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {getversion} from '../modules/getversionModules.nf' params(params)

// define workflow component
workflow vcfpredict {

    main:

      getversion()

    emit:

      software_json = gerversion.out.getversion_json

}

