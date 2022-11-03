#!/usr/bin/env python3

import json
import os
import sys

dirname = os.path.dirname(__file__)

scenarios = [
    ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'dryRun', 'NOW_ALIGN_TO_REF_dryRun', 'NOW_VARCALL_dryRun', 'CREATE_ANTIBIOGRAM_dryRun', 'yes', 'yes'],
    ['fastq', "*_R{1,2}.fastq.gz", 'null', 'null', 'fail', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'no', 'no'],
    ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'fail', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'no', 'no'],
    ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'dryRun', 'fail', 'null', 'null', 'null', 'null', 'null', 'null', 'no', 'no'],
    ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'dryRun', 'dryRun', 'null', 'null', 'null', 'null', 'null', 'null', 'no', 'no'],
    ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'null', 'null', 'null', 'null', 'null', 'no', 'no'],
    ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'dryRun', 'NOW_ALIGN_TO_REF_dryRun', 'NOW_VARCALL_dryRun', 'CREATE_ANTIBIOGRAM_dryRun', 'no', 'no'],
    ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'fail', 'null', 'null', 'null', 'no', 'no'],
    ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'dryRun', 'null', 'null', 'null', 'no', 'no'],
    ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'dryRun', 'NOW_ALIGN_TO_REF_dryRun', 'null', 'null', 'no', 'no'],
]


def main(output_dir, resources_dir):
    print(dirname)
    count = 0

    for scenario in scenarios:
        
        scenario_no = count + 1
        output_dir_scenario = os.path.join(output_dir, "scenario{}".format(scenario_no))
        report_dir_scenario = os.path.join(output_dir_scenario, "reports")
        file_name_scenario_inputs_json = "inputs.scenario{}.json".format(scenario_no)
        file_name_scenario_manifest_json = "MANIFEST.scenario{}.json".format(scenario_no)
        # write scenario inputs.json
        scenario_dict = {
                "filetype": scenario[0], 
                "pattern": scenario[1],
                "checkFqValidity_isok": scenario[2],
                "checkBamValidity_isok": scenario[3],
                "countReads_runfastp": scenario[4],
                "fastp_enoughreads": scenario[5],
                "kraken2_runmykrobe": scenario[6],
                "identifyBacContam_rundecontam": scenario[7],
                "downloadContamGenomes_fapass": scenario[8],
                "summary_doWeAlign": scenario[9],
                "alignToRef_doWeVarCall": scenario[10],
                "minos_isSampleTB": scenario[11], 
                "vcfmix": scenario[12],
                "gnomon": scenario[13],
                "output_dir": output_dir_scenario,
                "resources_dir": resources_dir,
                "report_dir": report_dir_scenario
                }
        
        with open(os.path.join(dirname, file_name_scenario_inputs_json), "w") as inputs_json:
            inputs_json.write(json.dumps(scenario_dict, indent=2))

        # write scenario MANIFEST.json
        manifest_dict = { 
                         "mainWorkflowURL": "main.nf",
                         "inputFileURLs": [
                              file_name_scenario_inputs_json
                              ],
                         "engineOptions": "-stub -config ./project/testing.config -profile agc"
                         }

        with open(os.path.join(dirname, file_name_scenario_manifest_json), "w") as manifest_json:
            manifest_json.write(json.dumps(manifest_dict, indent=2))
        
        count += 1

if __name__ == '__main__':
    output_dir = sys.argv[1]
    resources_dir = sys.argv[2]
    main(output_dir, resources_dir)
