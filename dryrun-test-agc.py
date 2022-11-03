#!/usr/bin/env python3

import json
import os
import sys
import subprocess
import argparse

parser = argparse.ArgumentParser(
    prog = 'dryrun-test-agc.py',
    description = 'For each test scenario generates inputs.json and MANIFEST.json files and submits a workflow run via Amazon Genomics CLI.'
)
parser.add_argument('-o', '--output_dir', dest='output_dir', help='parent output directory for test runs', required=True)
parser.add_argument('-r', '--resources_dir', dest='resources_dir', help='path to TB pipeline resources', required=True)
parser.add_argument('-c', '--context', dest='agc_context', help='AGC context to run the workflow in', default='ondemand', required=False)

scenarios = [
    ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'dryRun', 'NOW_ALIGN_TO_REF_dryRun', 'NOW_VARCALL_dryRun', 'CREATE_ANTIBIOGRAM_dryRun', 'yes', 'yes'],
    # ['fastq', "*_R{1,2}.fastq.gz", 'null', 'null', 'fail', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'no', 'no'],
    # ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'fail', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'no', 'no'],
    # ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'dryRun', 'fail', 'null', 'null', 'null', 'null', 'null', 'null', 'no', 'no'],
    # ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'dryRun', 'dryRun', 'null', 'null', 'null', 'null', 'null', 'null', 'no', 'no'],
    # ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'null', 'null', 'null', 'null', 'null', 'no', 'no'],
    # ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'dryRun', 'NOW_ALIGN_TO_REF_dryRun', 'NOW_VARCALL_dryRun', 'CREATE_ANTIBIOGRAM_dryRun', 'no', 'no'],
    # ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'fail', 'null', 'null', 'null', 'no', 'no'],
    # ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'dryRun', 'null', 'null', 'null', 'no', 'no'],
    # ['fastq', "*_R{1,2}.fastq.gz", 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'dryRun', 'NOW_ALIGN_TO_REF_dryRun', 'null', 'null', 'no', 'no'],
]

file_name_inputs_json = "inputs.json"
file_name_manifest_json = "MANIFEST.json"

def main(output_dir, resources_dir, context):

    print("submitting test workflow runs via Amazon Genomics CLI...")
    print("AGC context          : {}".format(context))
    print("output directory     : {}".format(output_dir))
    print("TB pipeline resources: {}".format(resources_dir))
    print("")
    print("backing up {} and {}...".format(file_name_inputs_json, file_name_manifest_json))
    
    # back up inputs.json and MANIFEST.json
    backup_files([file_name_inputs_json, file_name_manifest_json])

    # for each scenario...
    count = 0
    for scenario in scenarios:
        
        scenario_no = count + 1
        
        output_dir_scenario = os.path.join(output_dir, "scenario{}".format(scenario_no))
        report_dir_scenario = os.path.join(output_dir_scenario, "reports")
        scenario.append(output_dir_scenario)
        scenario.append(resources_dir)
        scenario.append(report_dir_scenario)
        
        # write scenario inputs.json
        write_inputs_json(scenario, file_name_inputs_json)

        # write scenario MANIFEST.json
        write_manifest_json(file_name_inputs_json, file_name_manifest_json)
        
        # submit workflow run
        print("submitting workflow run for scenario {}...".format(scenario_no))
        submit_workflow_run(context)
        
        count += 1
        break
    
    # restore inputs.json and MANIFEST.json
    print("restoring {} and {}...".format(file_name_inputs_json, file_name_manifest_json)) 
    restore_files([file_name_inputs_json, file_name_manifest_json])


def backup_files(file_names):

    for file_name in file_names:
        output = subprocess.check_output(
            "mv -v {} {}.bk".format(file_name, file_name),
            stderr=subprocess.STDOUT,
            shell=True
        )
        print(output.decode("utf-8").strip())
    

def restore_files(file_names):

    for file_name in file_names:
        output = subprocess.check_output(
            "mv -v {}.bk {}".format(file_name, file_name),
            stderr=subprocess.STDOUT,
            shell=True
        )
        print(output.decode("utf-8").strip())
    

def write_inputs_json(scenario, file_name):
    
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
        "output_dir": scenario[14],
        "resources_dir": scenario[15],
        "report_dir": scenario[16]
    }
        
    with open(file_name, "w") as inputs_json:
        inputs_json.write(json.dumps(scenario_dict, indent=2))    


def write_manifest_json(file_name_inputs_json, file_name):
    
    manifest_dict = {
        "mainWorkflowURL": "main.nf",
        "inputFileURLs": [
            file_name_inputs_json
        ],
        "engineOptions": "-stub -config ./project/testing.config -profile agc"
    }

    with open(file_name, "w") as manifest_json:
        manifest_json.write(json.dumps(manifest_dict, indent=2))
        

def submit_workflow_run(context):
    output = subprocess.check_output(
        "agc workflow run -c {} -n tbpipeline".format(context),
        stderr=subprocess.STDOUT,
        shell=True
    )
    print(output.decode("utf-8").strip())

    
if __name__ == '__main__':
    
    args = parser.parse_args()
    main(args.output_dir, args.resources_dir, args.agc_context)
