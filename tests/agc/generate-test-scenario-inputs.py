import json
import os
dirname = os.path.dirname(__file__)

count = 1
scenario1 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun',
                     'dryRun', 'NOW_ALIGN_TO_REF_dryRun', 'NOW_VARCALL_dryRun', 'CREATE_ANTIBIOGRAM_dryRun', 'yes', 'yes']
scenario2 = ['fastq', '"*_R{1,2}.fastq.gz"', 'null', 'null', 'fail', 'null', 'null', 'null', 'null', 'null', 'null',
                     'null', 'no', 'no']
scenario3 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'fail', 'null', 'null', 'null', 'null', 'null', 'null',
                     'null', 'no', 'no']
scenario4 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'fail', 'null', 'null', 'null', 'null', 'null',
                     'null', 'no', 'no']
scenario5 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'null', 'null', 'null', 'null', 'null',
                     'null', 'no', 'no']
scenario6 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'null', 'null', 'null', 'null',
                     'null', 'no', 'no']
scenario7 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun',
                     'dryRun', 'NOW_ALIGN_TO_REF_dryRun', 'NOW_VARCALL_dryRun', 'CREATE_ANTIBIOGRAM_dryRun', 'no', 'no']
scenario8 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun',
                     'fail', 'null', 'null', 'null', 'no', 'no']
scenario9 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun',
                     'dryRun', 'null', 'null', 'null', 'no', 'no']
scenario10 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun',
                      'dryRun', 'NOW_ALIGN_TO_REF_dryRun', 'null', 'null', 'no', 'no']

scenarios = []
for num in range(1, 11):
    scenario = 'scenario' + str(num)
    scenarios.append(locals()[scenario])

for scenario in scenarios:
    scenario_dict = {"filetype": scenario[0], "pattern": scenario[1],"checkFqValidity_isok": scenario[2],
     "checkBamValidity_isok": scenario[3],"countReads_runfastp": scenario[4],"fastp_enoughreads": scenario[5],
     "kraken2_runmykrobe":scenario[6],"identifyBacContam_rundecontam":scenario[7],
     "downloadContamGenomes_fapass":scenario[8],"summary_doWeAlign":scenario[9],"alignToRef_doWeVarCall": scenario[10],
     "minos_isSampleTB":scenario[11],"vcfmix":scenario[12],
     "gnomon":scenario[13],"output_dir":"s3://agc-711700981500-eu-west-2/project/tbpipeline/tests/","resources_dir": "s3://tbpipeline/resources"}
    with open(os.path.join(dirname, "output/" + "inputs.scenario" + str(count) + ".json"), "w") as file_handle:
        file_handle.write(json.dumps(scenario_dict))

    manifest_dict = {"mainWorkflowURL": "main.nf","inputFileURLs": [
    "inputs.json"],"engineOptions": "-stub -config ./project/testing.config -profile agc"}

    with open(os.path.join(dirname, "output/" + "MANIFEST" + str(count) + ".json"), "w") as file_handle:
        file_handle.write(json.dumps(manifest_dict))
    count += 1
