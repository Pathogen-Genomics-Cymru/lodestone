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
      input_seqs_json

    main:

      alignToRef(input_seqs_json)

      callVarsMpileup(alignToRef.out.alignToRef_bam)

      callVarsCortex(alignToRef.out.alignToRef_bam)

      minos(alignToRef.out.alignToRef_bam.join(callVarsCortex.out.cortex_vcf, by: 0).join(callVarsMpileup.out.mpileup_vcf, by: 0))

      gvcf(alignToRef.out.alignToRef_bam.join(minos.out.minos_vcf, by: 0))
}
