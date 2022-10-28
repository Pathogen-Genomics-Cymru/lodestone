import json
import os
import sys
import argparse

# function to read kraken report
def read_kraken_report(input, pct_threshold, num_threshold):
    # input - full path to kraken report
    # pct_threshold - min coverage, as %%
    # num_threshold - min coverage, as no. of reads

    # the Kraken report is assumed to be the standard format: 6 tab-delimited columns, with one line per taxon. This is described at https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats. We will confirm this as we parse.
    S = []
    G = []
    G1 = []
    F = []
    non_human_species_detected = 0

    # open file and read it line by line
    lineCount = 0
    with open(input, "r") as f:
        for l in f:
            line = l
            lineCount += 1
            # trim and split line based on \t
            line = line.strip().split("\t")
            # skip if not 6 columns. In a correctly formatted Kraken report, each row will have 6 columns; we will skip those that don't
            if len(line) != 6:
                continue
            pc_frags = line[0] # defined as "percentage of fragments covered by the clade rooted at this taxon. NOTE: for this purpose, 'fragment' is synonymous with 'read'.
            num_frags_rooted = line[1] # defined as "number of fragments covered by the clade rooted at this taxon", i.e. all fragments at this taxonomic level AND LOWER
            num_frags_direct = line[2] # defined as "number of fragments assigned directly to this taxon"
            rank_code = line[3] # defined as "rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies". Takes the form of one letter, optionally followed by one number.
            ncbi_taxon_id = line[4] # defined as "NCBI taxonomic ID number"
            name = line[5] # defined as "scientific name"

            # convert strings to float or int
            pc_frags = float(pc_frags.strip())
            num_frags_rooted = int(num_frags_rooted.strip())
            num_frags_direct = int(num_frags_direct.strip())
            ncbi_taxon_id = ncbi_taxon_id.strip()
            rank_code = rank_code.strip()
            name = name.strip()

            # skip classifications not supported by a min % of fragments
            if (pc_frags < pct_threshold) & (name != 'Homo sapiens'):
                continue
            # skip classifications not supported by a min no. of fragments
            if (num_frags_rooted < num_threshold) & (name != 'Homo sapiens'):
                continue

            if (isinstance(pc_frags, float)) & (isinstance(num_frags_rooted, int)) & (isinstance(num_frags_direct, int)) & (isinstance(rank_code, str)) & (isinstance(ncbi_taxon_id, str)) & (isinstance(name, str)):
                if rank_code == 'S':
                    S.append([num_frags_rooted, pc_frags, name, ncbi_taxon_id])
                    if name != 'Homo sapiens': non_human_species_detected += 1
                elif rank_code == 'G':
                    G.append([num_frags_rooted, pc_frags, name, ncbi_taxon_id])
                elif rank_code == 'F':
                    F.append([num_frags_rooted, pc_frags, name, ncbi_taxon_id])

                # Kraken does not resolve classifications among the Mycobacteriaceae as well as Mykrobe. At best, it can detect species complexes. We shall retain these classifications to look at later, as they may indicate whether this is a mixed-mycobacterial sample.
                if (name.startswith('Mycobact')) & (rank_code == 'G1'):
                    G1.append([num_frags_rooted, pc_frags, name, ncbi_taxon_id])
            else:
                sys.exit('ERROR: malformatted Kraken report, at line %d' %(lineCount))
    
    return S, G, G1, F, non_human_species_detected

# define output function
def parse_kraken_report(S, G, G1, F, non_human_species_detected, pct_threshold, num_threshold):
    # arguments are the output from read_kraken_report function

    # define warnings lists
    warnings = []

    # define output
    out = {
        "Thresholds": {
            "reads": num_threshold,
            "percentage": pct_threshold
        }
    }

    if (len(S) == 0) | (non_human_species_detected == 0):
        warnings.append("warning: no species classifications meet thresholds of > %d reads and > %s %% of total reads (human excepted)" %(num_threshold, pct_threshold))
    if len(G) == 0:
        warnings.append("warning: no genus classifications meet thresholds of > %d reads and > %s %% of total reads" %(num_threshold, pct_threshold))
    if len(F) == 0:
        warnings.append("warning: no family classifications meet thresholds of > %d reads and > %s %% of total reads" %(num_threshold, pct_threshold))

    top_family = ''
    no_of_reads_assigned_to_top_family = ''
    top_genus = ''
    top_species = ''
    contaminant_species_found = 0
    contaminant_mycobacterium_found = 0

    for x in range (0, 3):
        arr = []
        clade = ''
        if x == 0:
            arr = F
            clade = 'Family'
        elif x == 1:
            arr = G
            clade = 'Genus'
        elif x == 2:
            arr = S
            clade = 'Species'

        if (len(arr) == 0):
            continue
        
        sorted_arr = arr.copy()
        sorted_arr.sort(reverse=True)
        if x == 0:
            top_family = sorted_arr[0][2]
            no_of_reads_assigned_to_top_family = sorted_arr[0][0]
        elif x == 1:
            top_genus = sorted_arr[0][2]
        elif x == 2:
            top_species = sorted_arr[0][2]
        
        for y in range (0, len(sorted_arr)):
            hash = {
                'reads': sorted_arr[y][0],
                'percentage': sorted_arr[y][1],
                'name': sorted_arr[y][2],
                'taxon_id': sorted_arr[y][3]
            }
            if clade not in out: out[clade] = []
            out[clade].append(hash)
            if (x == 2) & (y > 0) & (sorted_arr[y][2] != 'Homo sapiens'):
                contaminant_species_found += 1 # raise a warning if a non-human species is detected that is NOT the top hit, as this indicates the sample is mixed or contaminated
                if (sorted_arr[y][2].startswith('Mycobact')):
                    contaminant_mycobacterium_found += 1
    
    if contaminant_species_found > 0:
        if contaminant_mycobacterium_found > 0:
            warnings.append("warning: sample is mixed or contaminated (contains reads from multiple non-human species). Contaminants (i.e. minority species) include one or more mycobacteria. Defer to Mykrobe report for superior mycobacterial classification")
        else:
            warnings.append("warning: sample is mixed or contaminated (contains reads from multiple non-human species)")
    if (top_family.startswith("Mycobact")) & ((top_genus.startswith("Mycobact") == False) | (top_species.startswith("Mycobact") == False)):
        warnings.append("warning: top family classification is mycobacterial, but this is not consistent with top genus and species classifications")

    # IF THE TOP FAMILY IS MYCOBACTERIACEAE (WHICH CAN ONLY BE THE CASE IF MINIMUM COVERAGE THRESHOLDS ARE MET), WE WILL ALSO REPORT THE KRAKEN 'G1' CLASSIFICATIONS. THESE MAY INDICATE WHETHER THIS IS A MIXED MYCOBACTERIAL SAMPLE.
    if top_family == "Mycobacteriaceae":
        if no_of_reads_assigned_to_top_family < 100000:
            if "Errors" not in out: out['Errors'] = []
            out['Errors'].append("error: there are < 100k reads classified as Mycobacteriaceae")
            out['Mykrobe'] = 'false'
        else:
            out['Mykrobe'] = 'true' # as the sample is predominantly Mycobacteriaceae, we recommend the user invoke Mykrobe for higher-resolution classification. Later in the workflow, we will be using this value in a text comparison. It MUST be lower-case here otherwise it will be mistaken for a Boolean (TRUE/FALSE) instead.

        if (len(G1) == 0):
            warnings.append("warning: top family is Mycobacteriaceae but no G1 (species complex) classifications meet thresholds of > %d reads and > %s %% of total reads (this is not necessarily a concern as not all mycobacteria have this taxonomy)" %(num_threshold, pct_threshold))
        else:
            sorted_G1 = G1.copy()
            sorted_G1.sort(reverse=True)
            if "Species complex" not in out: out['Species complex'] = []
            for x in range(0, len(sorted_G1)):
                hash = {
                    'reads': sorted_G1[x][0],
                    'percentage': sorted_G1[x][1],
                    'name': sorted_G1[x][2],
                    'taxon_id': sorted_G1[x][3]
                }
                out['Species complex'].append(hash)
            
            if len(sorted_G1) > 1:
                warnings.append("warning: sample contains multiple mycobacterial species complexes (for superior classification of mixed mycobacteria, defer to Mykrobe report)")
    else:
        if "Errors" not in out: out['Errors'] = []
        out['Errors'].append("error: top family is not Mycobacteriaceae")
        out['Mykrobe'] = 'false'

    if len(warnings) == 0:
        warnings.append('')
    out['Warnings'] = warnings

    return out

# process requirements function (check input)
def process_requirements(args):
    # REQUIREMENTS
    in_file = args[1]
    out_file = args[2]
    pct_threshold = float(args[3])
    num_threshold = int(args[4])

    # check if input file exists
    if not os.path.exists(in_file):
        sys.exit('ERROR: cannot find %s' %(in_file))

    # check if input file is empty
    if os.stat(in_file).st_size == 0:
        sys.exit('ERROR: %s is empty' %(in_file))

    # check if output file ends with .json
    if not out_file.endswith('.json'):
        sys.exit('ERROR: output file %s must end with suffix .json' %out_file)
    
    # check if pct_threshold and num_threshold are positive floar/int
    if not (isinstance(num_threshold, int) & (num_threshold > 0)):
        sys.exit('ERROR: %d is not a positive integer' %(num_threshold))
    if not (isinstance(pct_threshold, float) & (pct_threshold > 0)):
        sys.exit('ERROR: %f is not a positive number' %(pct_threshold))

    # check if pct_threshold is more than 100
    if pct_threshold > 100:
        sys.exit('ERROR: %f is a %% and cannot be > 100' %(pct_threshold))

    return

# call main function
if __name__ == "__main__":

    # set command line arguments
    description = 'This script will parse a Kraken output file to report all family/genus/species classifications in the sample, plus species complex classifications if the dominant family is Mycobacteriaceae.\n'
    description += 'We require a min. coverage of x%% of the total reads AND a min. number of reads PER CLASSIFICATION. Set these to 0 to report everything.\n'
    usage = 'python parse_kraken_report2.py [path to Kraken report] [path to output file; must end .json] [min. coverage, as %%] [min. coverage, as no. of reads]\n'
    usage += 'E.G.:\tpython parse_kraken_report2.py report.txt out.json 1 10000\n\n\n'
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('in_file', metavar='in_file', type=str, help='Path to Kraken report')
    parser.add_argument('out_file', metavar='out_file', type=str, help='Path to output file; must end .json')
    parser.add_argument('pct_threshold', metavar='pct_threshold', type=float, help='Min. coverage, as %%')
    parser.add_argument('num_threshold', metavar='num_threshold', type=int, help='Min. coverage, as no. of reads. Should be a positive integer.')

    args = parser.parse_args()

    # requirements
    process_requirements(sys.argv)
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    pct_threshold = float(sys.argv[3])
    num_threshold = int(sys.argv[4])

    # read kraken report
    S = []
    G = []
    G1 = []
    F = []
    non_human_species_detected = 0
    S, G, G1, F, non_human_species_detected = read_kraken_report(in_file, pct_threshold, num_threshold)

    # parse kraken report and generate output
    out = parse_kraken_report(S, G, G1, F, non_human_species_detected, pct_threshold, num_threshold)
    
    # CREATE OUTPUT FILE
    with open(out_file, 'w') as f:
        json.dump(out, f, indent = 4)