// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {alignToRef} from '../modules/clockworkModules.nf' params(params)
include {callVarsMpileup} from '../modules/clockworkModules.nf' params(params)

// define workflow component
workflow clockwork {

    take:
      input_seqs
      json

    main:

      alignToRef(input_seqs, json)
      callVarsMpileup(alignToRef.out.alignToRef_mpileup) 
}
