// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {alignToRef} from '../modules/clockworkModules.nf' params(params)
include {callVarsMpileup} from '../modules/clockworkModules.nf' params(params)
include {callVarsCortex} from '../modules/clockworkModules.nf' params(params)
include {minos} from '../modules/clockworkModules.nf' params(params)
include {gvcf} from '../modules/clockworkModules.nf' params(params)
include {getRefFromJSON} from '../modules/clockworkModules.nf' params(params)
include {getRefCortex} from '../modules/clockworkModules.nf' params(params)
         
// define workflow component
workflow clockwork {

    take:
      input_seqs_json

    main:
      //get just the json
      json = input_seqs_json.map{it[4]}
      do_we_align = input_seqs_json.map{it[5]}
      sample_name = input_seqs_json.map{it[0]}
      
      getRefFromJSON(json, do_we_align, sample_name)
      alignToRef(input_seqs_json, getRefFromJSON.out)
      

      callVarsMpileup(alignToRef.out.alignToRef_bam)

      getRefCortex(alignToRef.out.alignToRef_bam)
      callVarsCortex(alignToRef.out.alignToRef_bam, getRefCortex.out)

      minos(alignToRef.out.alignToRef_bam.join(callVarsCortex.out.cortex_vcf, by: 0).join(callVarsMpileup.out.mpileup_vcf, by: 0))

      gvcf(alignToRef.out.alignToRef_bam.join(minos.out.minos_vcf, by: 0))

    emit:

      mpileup_vcf = callVarsMpileup.out.mpileup_vcf.join(minos.out.minos_report, by: 0)
      minos_vcf = minos.out.minos_vcf.join(alignToRef.out.alignToRef_report, by: 0)
      reference = getRefFromJSON.out
      bam = alignToRef.out.alignToRef_bam

}
