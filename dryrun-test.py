#!/usr/bin/env python3

import subprocess
import filecmp

def run_tests():

    scenario1 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'dryRun', 'NOW_ALIGN_TO_REF_dryRun', 'NOW_VARCALL_dryRun']
    scenario2 = ['fastq', '"*_R{1,2}.fastq.gz"', 'null', 'null', 'fail', 'null', 'null', 'null', 'null', 'null', 'null']
    scenario3 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'fail', 'null', 'null', 'null', 'null', 'null', 'null']
    scenario4 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'fail', 'null', 'null', 'null', 'null', 'null']
    scenario5 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'null', 'null', 'null', 'null', 'null']
    scenario6 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'null', 'null', 'null', 'null']
    scenario7 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'dryRun', 'NOW_ALIGN_TO_REF_dryRun', 'NOW_VARCALL_dryRun']
    scenario8 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'fail', 'null', 'null']
    scenario9 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'dryRun', 'null', 'null']
    scenario10 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'dryRun', 'NOW_ALIGN_TO_REF_dryRun', 'null']
    
    scenarios=[]
    for num in range(1, 11):
        scenario = 'scenario' + str(num)
        scenarios.append(locals()[scenario])

    count=1
    for scenario in scenarios:
        toRun = 'NXF_VER=20.11.0-edge nextflow run main.nf -stub -config testing.config' + ' --filetype ' + scenario[0] + ' --pattern ' + scenario[1] \
        + ' --checkFqValidity_isok ' + scenario[2] + ' --checkBamValidity_isok ' + scenario[3] + ' --countReads_runfastp ' + scenario[4] \
        + ' --fastp_enoughreads ' + scenario[5] + ' --kraken2_runmykrobe ' + scenario[6] + ' --identifyBacContam_rundecontam ' + scenario[7] \
        + ' --downloadContamGenomes_fapass ' + scenario[8] + ' --summary_doWeAlign ' + scenario[9] + ' --alignToRef_doWeVarCall ' + scenario[10] \
        + ' > scenario' + str(count) + '.txt' 
        result = subprocess.run([toRun], shell=True)
        count+=1

def compare_with_truth_set():

    for num in range(1, 11):
        filename = 'scenario' + str(num) + '.txt'
        truthset = './tests/' + filename 
        filecmp.cmp(filename, truthset)
        if filecmp.cmp:
            print (filename + ' pass')
        else:
            print (filename + ' fail')

def main():

    run_tests()
    compare_with_truth_set()

if __name__ == '__main__':
    main()

