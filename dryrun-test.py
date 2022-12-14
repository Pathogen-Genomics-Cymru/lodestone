#!/usr/bin/env python3

import subprocess
import os

file_dir = os.path.dirname(os.path.abspath(__file__))
truth_dir = os.path.join(file_dir, 'tests')

def go():

    scenario1 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'dryRun', 'NOW_ALIGN_TO_REF_dryRun', 'NOW_VARCALL_dryRun', 'CREATE_ANTIBIOGRAM_dryRun', 'yes', 'yes']
    scenario2 = ['fastq', '"*_R{1,2}.fastq.gz"', 'null', 'null', 'fail', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'no', 'no']
    scenario3 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'fail', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'no', 'no']
    scenario4 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'fail', 'null', 'null', 'null', 'null', 'null', 'null', 'no', 'no']
    scenario5 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'null', 'null', 'null', 'null', 'null', 'null', 'no', 'no']
    scenario6 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'null', 'null', 'null', 'null', 'null', 'no', 'no']
    scenario7 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'dryRun', 'NOW_ALIGN_TO_REF_dryRun', 'NOW_VARCALL_dryRun', 'CREATE_ANTIBIOGRAM_dryRun', 'no', 'no']
    scenario8 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'fail', 'null', 'null', 'null', 'no', 'no']
    scenario9 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'dryRun', 'null', 'null', 'null', 'no', 'no']
    scenario10 = ['fastq', '"*_R{1,2}.fastq.gz"', 'OK', 'null', 'dryRun', 'dryRun', 'dryRun', 'NOW_DECONTAMINATE_dryRun', 'dryRun', 'NOW_ALIGN_TO_REF_dryRun', 'null', 'null', 'no', 'no']
    
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
        + ' --minos_isSampleTB ' + scenario[11] + ' --vcfmix ' + scenario[12] + ' --gnomonicus ' + scenario[13] \
        + ' > scenario' + str(count) + '.txt' 
        result = subprocess.run([toRun], shell=True)
        count+=1

    for num in range(1, 11):
        filename = 'scenario' + str(num) + '.txt'
        with open(filename) as in_file:
            data = in_file.readlines()
            tail = data[-28:]
        with open(filename, 'w') as out_file:
            out_file.write(''.join(tail))
        filesize = os.stat(filename).st_size
        truthset = os.path.join(truth_dir, filename)
        truthsize = os.stat(truthset).st_size
        if filesize == truthsize:
            print ('scenario' + str(num) + ' passed dry run')
        else:
            raise ValueError('scenario' + str(num) + ' failed dry run')

def main():
    go()

if __name__ == '__main__':
    main()
