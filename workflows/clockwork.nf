// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {alignToRef} from '../modules/clockworkModules.nf' params(params)
include {callVarsMpileup} from '../modules/clockworkModules.nf' params(params)
include {callVarsCortex} from '../modules/clockworkModules.nf' params(params)
include {minos} from '../modules/clockworkModules.nf' params(params)
include {gvcf} from '../modules/clockworkModules.nf' params(params)

// define workflow component
workflow clockwork {

    take:
      input_seqs
      json

    main:

      alignToRef(input_seqs, json)
      callVarsMpileup(alignToRef.out.alignToRef_bam)
      callVarsCortex(alignToRef.out.alignToRef_bam)
      minos(alignToRef.out.alignToRef_bam, callVarsCortex.out.cortex_vcf, callVarsMpileup.out.mpileup_vcf)
      gvcf(alignToRef.out.alignToRef_bam, minos.out.minos_vcf)
}
