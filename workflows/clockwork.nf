// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {alignToRef} from '../modules/clockworkModules.nf'
include {callVarsMpileup} from '../modules/clockworkModules.nf'
include {callVarsCortex} from '../modules/clockworkModules.nf'
include {minos} from '../modules/clockworkModules.nf'
include {gvcf} from '../modules/clockworkModules.nf'
include {getRefFromJSON} from '../modules/clockworkModules.nf'
include {getRefCortex} from '../modules/clockworkModules.nf'
         
// define workflow component
workflow clockwork {

    take:
      input_seqs_json //tuple of sample_name, fq1, fq2, software_version_json, 
                      //speciation_json, do_we_align?

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

      minos(alignToRef.out.alignToRef_bam
            .join(callVarsCortex.out.cortex_vcf, by: 0)
            .join(callVarsMpileup.out.mpileup_vcf, by: 0))

      gvcf(alignToRef.out.alignToRef_bam.join(minos.out.minos_vcf, by: 0))

      report_for_ntm = gvcf.out.gvcf_report_resistance
      sample_and_fqs = input_seqs_json.map{it[0,1,2]}
      profiler_input_fq = sample_and_fqs.join(report_for_ntm, by:0)

    emit:
      profiler_input_vcf = gvcf.out.tbprofiler
      profiler_input_fq = profiler_input_fq
}
