#!/usr/bin/env python3
 import json
import os
import sys
import argparse
import re
import copy

# process requirements function
def process_requirements(args):
    stats_file = args[1]
    report_file = args[2]

    # check if input files exist
    if not os.path.exists(stats_file):
        sys.exit('ERROR: cannot find %s' %(stats_file))

    if not os.path.exists(report_file):
        sys.exit('ERROR: cannot find %s' %(report_file))

    # check if input files are empty
    if os.stat(stats_file).st_size == 0:
        sys.exit('ERROR: %s is empty' %(stats_file))

    if os.stat(report_file).st_size == 0:
        sys.exit('ERROR: %s is empty' %(report_file))

    # check IDs from the file names
    res_stats1 = re.findall(r"^.+[\/|\\](.*?)\_alignmentStats\.json$", stats_file)
    res_stats2 = re.findall(r"^(.*?)\_alignmentStats\.json$", stats_file)
    if len(res_stats1) > 0:
        sample_id_STA = res_stats1[0]
    elif len(res_stats2) > 0:
        sample_id_STA = res_stats2[0]
    else:
        sample_id_STA = ''

    res_report1 = re.findall(r"^.+[\/|\\](.*?)\_species\_in\_sample\_pass\_(.*)\.json$", report_file)
    res_report2 =  re.findall(r"^(.*?)_species\_in\_sample\_pass\_(.*)\.json$", report_file)
    res_report3 = re.findall(r"^.+[\/|\\](.*?)\_species\_in\_sample\.json$", report_file)
    res_report4 = re.findall(r"^(.*?)_species\_in\_sample\.json$", report_file)

    if len(res_report1) > 0:
        sample_id_REP = res_report1[0][0]
    elif len(res_report2) > 0:
        sample_id_REP = res_report2[0][0]
    elif len(res_report3) > 0:
        sample_id_REP = res_report3[0]
    elif len(res_report4) > 0:
        sample_id_REP = res_report4[0]
    else:
        sample_id_REP = ''
    print(sample_id_REP)

    sample_id = ''
    if sample_id_STA != sample_id_REP:
        sys.exit("ERROR: the sample IDs of %s and %s are mismatched" %(stats_file, report_file))
    else:
        sample_id = sample_id_STA

    if sample_id == '':
        sys.exit("ERROR: could not identify sample ID from the filename of either %s or %s" %(stats_file, report_file))

    return sample_id

# function to read input files
def read_and_parse_input_files(stats_file, report_file):
    # READ INPUT FILES
    with open(stats_file, 'r') as f:
        stats = json.load(f) 
    with open(report_file, 'r') as f:
        report = json.load(f)
    
    # WHAT % OF THE REFERENCE GENOME IS COVERED AT 10-FOLD DEPTH?
    pc_10fold = stats['bam_stats']['coverage_breadth']['10']
    # HOW MANY READS WERE MAPPED TO THE REFERENCE GENOME?
    num_mapped = stats['bam_stats']['reads_mapped']
    # WHAT IS THE AVERAGE READ MAPPING QUALITY?
    avg_mapq = stats['bam_stats']['average_quality']

    # ARE THERE ANY ERRORS?
    errors = []
    num_errors = 0
    if num_mapped < 100000:
        num_errors += 1
        errors.append("error: < 100k reads could be mapped to the reference genome")
    if pc_10fold < 50:
        num_errors += 1
        errors.append("error: < 50%% of the reference genome is covered at 10-fold depth")
    if avg_mapq < 10:
        num_errors += 1
        errors.append("error: alignments to the reference genome have average mapping quality < 10")

    # IF THERE ARE NO EXISTING WARNINGS BUT WE'VE NOW RULED THERE ARE ERRORS, WE REVISE THE WARNING MESSAGE TO POINT THIS OUT
    warnings = report['warnings']
    num_warnings = len(warnings)
    no_warning_message = 0
    for warning in warnings:
        if (warning == 'no warnings raised'):
            no_warning_message += 1
    
    if ((no_warning_message == 1) & (num_warnings == 1) & (num_errors > 0)):
        warnings = []
        if num_errors == 1:
            warnings.append("there was %d error but no warnings" %num_errors)
        elif num_errors > 1:
            warnings.append("there was %d errors but no warnings" %num_errors)
        report['warnings'] = warnings
    
    out = copy.deepcopy(report)
    bam_stats = stats['bam_stats']
    out['bam_stats'] = bam_stats
    if num_errors > 0:
        out['errors'] = errors
        out['summary_questions']['continue_to_clockwork'] = 'no'

    return out

# call main function
if __name__ == "__main__":
    # set command line arguments
    description = 'This script will produce the test_report.json\n'
    usage = 'python create_final_json.py [path to stats json] [path to report json]'
    usage += 'E.G.:\tpython create_final_json.py stats.json report.json\n\n\n'
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('stats_json', metavar='stats_json', type=str, help='Path to stats json')
    parser.add_argument('report_json', metavar='report_json', type=str, help='Path to report json')

    args = parser.parse_args()

    # requirements
    sample_id = process_requirements(sys.argv)
    stats_file = sys.argv[1]
    report_file = sys.argv[2]

    # read and parse input files
    out = read_and_parse_input_files(stats_file, report_file)

    # generate output file
    output_path = sample_id + "_report.json"
    with open(output_path, "w") as f:
        json.dump(out, f, indent = 4)

