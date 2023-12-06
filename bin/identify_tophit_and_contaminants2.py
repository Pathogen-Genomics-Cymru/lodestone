#!/usr/bin/env python3

import json
import os
import sys
import argparse
import re
import copy

# define process requirements function
def process_requirements(args):
    # REQUIREMENTS
    afanc_json = args[1]
    kraken_json = args[2]
    assembly_file = args[3]
    supposed_species = args[4]
    unmix_myco = args[5]
    myco_dir = args[6]
    prev_species_json = args[7]
    
    """
    # check if input files exist and not empty
    if not os.path.exists(afanc_json):
        sys.exit('ERROR: cannot find %s' %(afanc_json))
    if os.stat(afanc_json).st_size == 0:
        sys.exit('ERROR: %s is empty' %(afanc_json))

    if not os.path.exists(kraken_json):
        sys.exit('ERROR: cannot find %s' %(kraken_json))
    if os.stat(kraken_json).st_size == 0:
        sys.exit('ERROR: %s is empty' %(kraken_json))

    if not os.path.exists(assembly_file):
        sys.exit('ERROR: cannot find %s' %(assembly_file)) # from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
    if os.stat(assembly_file).st_size == 0:
        sys.exit('ERROR: %s is empty' %(assembly_file))

    if not os.path.exists(myco_dir) and not bucket_exists(myco_dir):
        sys.exit('ERROR: cannot find %s' %(myco_dir))

    if (prev_species_json != 'null'):
        if not os.path.exists(prev_species_json):
            sys.exit('ERROR: cannot find %s' %(prev_species_json))
        if os.stat(prev_species_json).st_size == 0:
            sys.exit('ERROR: %s is empty' %(prev_species_json))
     """
    
    species = ['abscessus', 'africanum', 'avium', 'bovis', 'chelonae', 'chimaera', 'fortuitum', 'intracellulare', 'kansasii', 'tuberculosis']
    for spec in species:
        spec_fasta_path = os.path.join(myco_dir, spec + '.fasta')
        spec_mmi_path = os.path.join(myco_dir, spec + '.mmi')

        """
        if myco_dir.startswith("s3://"):
            s3_myco_dir = myco_dir.replace("s3://", "")
            spec_fasta = s3_myco_dir.split("/", 1)[-1] + "/" + spec + ".fasta"
            s3_myco_dir =  s3_myco_dir.split("/", 1)[0]

            if not is_file_in_s3(s3_myco_dir, spec_fasta):
                sys.exit('ERROR: cannot find %s' %(spec_fasta_path))
        else:
            if not os.path.exists(spec_fasta_path):
                sys.exit('ERROR: cannot find %s' %(spec_fasta_path))

        if myco_dir.startswith("s3://"):
            s3_myco_dir = myco_dir.replace("s3://", "")
            spec_mmi = s3_myco_dir.split("/", 1)[-1] + "/" + spec + ".mmi"
            s3_myco_dir =  s3_myco_dir.split("/", 1)[0]

            if not is_file_in_s3(s3_myco_dir, spec_mmi):
                sys.exit('ERROR: cannot find %s' %(spec_mmi_path))
        else:
            if not os.path.exists(spec_fasta_path):
                sys.exit('ERROR: cannot find %s' %(spec_mmi_path))
        """
        
    if ((supposed_species != 'null') & (supposed_species not in species)):
        sys.exit('ERROR: if you provide a species ID, it must be one of either: abscessus|africanum|avium|bovis|chelonae|chimaera|fortuitum|intracellulare|kansasii|tuberculosis')

    if ((unmix_myco != 'yes') & (unmix_myco != 'no')):
        sys.exit('ERROR: \'unmix myco\' should be either \'yes\' or \'no\'')

    ## check IDs from the file names

    # get ID of a Mykrobe/Afanc report depending on which report is provided
    if afanc_json.endswith("_mykrobe_report.json"):
        sample_id_MYK = os.path.basename(afanc_json).split("_mykrobe")[0]
    elif afanc_json.endswith("_afanc_report.json"):
        sample_id_MYK = os.path.basename(afanc_json).split("_afanc")[0]
    else:
        sample_id_MYK = ''

    # get ID of a Kraken report
    if kraken_json.endswith("_kraken_report.json"):
        sample_id_KRA = os.path.basename(kraken_json).split("_kraken")[0]
    else:
        sample_id_KRA = ''

    # check if Mykrobe/Afanc report ID matches Kraken report ID
    sample_id = ''
    if sample_id_MYK != sample_id_KRA:
        sys.exit("ERROR: the sample IDs of %s and %s are mismatched" %(afanc_json, kraken_json))
    else:
        sample_id = sample_id_MYK

    # if previous top-hit/contaminant report is provided, find its ID and check if that matches sample ID identified earlier
    if prev_species_json != 'null':
        if prev_species_json.endswith("_species_in_sample_previous.json"):
            sample_id_PRE = os.path.basename(prev_species_json).split("_species")[0]
        else:
            sample_id_PRE = ''

        if sample_id != sample_id_PRE:
            sys.exit("ERROR: sample ID of the previous species JSON (%s) does not match the sample ID we have from the Kraken and afanc reports (%s)" %(prev_species_json, sample_id))

    # if sample ID could not be identified, produce an error message
    if sample_id == '':
        sys.exit("ERROR: could not identify sample ID from the filename of either %s or %s" %(afanc_json, kraken_json))

    return sample_id

# PARSE THE 'ASSEMBLY SUMMARY' FILE TO STORE URLs FOR EACH GENOME ASSOCIATED WITH A GIVEN TAXON ID. WE REQUIRE THAT THIS IS THE LATEST VERSION OF THE GENOME AVAILABLE IN REFSEQ. WE WILL LATER CHECK THAT (A) IT IS DEFINED AS A "COMPLETE GENOME", AND (B) IT HAS FULL GENOME REPRESENTATION.
## RefSeq defines a 'complete genome' (see ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt) as follows: "all chromosomes are gapless and have no runs of 10 or more ambiguous bases (Ns), there are no unplaced or unlocalized scaffolds, and all the expected chromosomes are present (i.e. the assembly is not noted as having partial genome representation). Plasmids and organelles may or may not be included in the assembly but if present then the sequences are gapless"
## RefSeq defines 'genome representation' (also in the above readme) as "whether the goal for the assembly was to represent the whole genome or only part of it." This takes one of two values: 'full' or 'partial'
def read_assembly_summary(assembly_file_path):
    urls = {}
    tax_ids = {}
    # open file and read it line by line
    with open(assembly_file_path, "r", encoding="utf-8") as f:
        for l in f:
            line = l
            if line.startswith("#"): continue
            # trim and split line based on \t
            line = line.strip().split("\t")
            refseq_category = line[4]
            species_taxid = line[6]
            ## strips out . characters from the species name
            species_name = line[7].replace(".", "")
            infraspecific_name = line[8]
            version_status = line[10]
            assembly_level = line[11]
            genome_rep = line[13]
            ftp_path = line[19]

            species_taxid_regex_res = re.findall(r"\d+", species_taxid)
            if ((version_status != 'latest') | (len(species_taxid_regex_res) == 0) | (ftp_path == 'na')): continue

            species_taxid = int(species_taxid)
            tax_ids[species_name] = species_taxid
            if species_taxid not in urls: urls[species_taxid] = []
            urls[species_taxid].append([species_name, infraspecific_name, ftp_path, assembly_level, genome_rep, refseq_category])

    return urls, tax_ids

## Function to match species according to the NCBI taxonomy for Mycobacteriaceae. Check is case-independent.
# Taxonomy: "Mycobacterium, Mycobacteroides, Mycolicibacter, Mycolicibacterium, and Mycolicibacillus"
def match_taxonomy(spec):
    if spec.lower().startswith('mycobact') or spec.lower().startswith('mycolicibac'):
        return True
    else:
        return False

# define main function to process data
def process_reports(afanc_json_path, kraken_json_path, supposed_species, unmix_myco, myco_dir_path, prev_species_json_path, urls, tax_ids, sample_id):

    # DEFINE OUTPUT
    out = {}
    warnings = []

    # OPEN JSON FILES
    with open(afanc_json_path, 'r') as f:
        afanc = json.load(f)
    with open(kraken_json_path, 'r') as f:
        kraken = json.load(f)
    prev_species = ''
    if (prev_species_json_path != 'null'):
        with open(prev_species_json_path, 'r') as f:
            prev_species = json.load(f)

    # WHAT IS THE TOP HIT MYCOBACTERIAL SPECIES IN THE SAMPLE, ACCORDING TO AFANC, AND ON THE BASIS OF % COVERAGE?
    species = []
    afanc_finds_nothing = 0

    for spec in afanc[sample_id]['phylogenetics']['species']:
        pc_coverage = afanc[sample_id]['phylogenetics']['species'][spec]['percent_coverage']
        median_depth = afanc[sample_id]['phylogenetics']['species'][spec]['median_depth']
        spec = spec.replace("_", " ")
        species.append([pc_coverage, median_depth, spec])
        if spec == 'Unknown': afanc_finds_nothing += 1

    sorted_species = species.copy()
    sorted_species.sort(reverse=True)
    pc_coverage_of_top_species = sorted_species[0][0]
    depth_of_top_species = sorted_species[0][1]
    top_species = sorted_species[0][2]

    out['top_hit'] = {}
    out['top_hit']['name'] = top_species
    num_afanc_species = len(sorted_species)

    # WE ARE ASSUMING THE TOP AFANC HIT, ON THE BASIS OF % COVERAGE, IS *ALSO* THE TOP HIT ON THE BASIS OF MEDIAN DEPTH. LET'S CONFIRM THIS, AND WARN IF THIS IS NOT THE CASE.
    # IT MAY BE POSSIBLE THAT A SAMPLE IS A MIXTURE OF SPECIES X (99% COVERAGE AT 10-FOLD DEPTH) AND SPECIES Y (98% COVERAGE AT 11-FOLD DEPTH). IN THIS CASE, ON WHAT BASIS DO WE CHOOSE A TOP HIT, GIVEN WE HAVE TO CHOOSE *ONE*?
    for x in range(1, len(sorted_species)):
        pc_coverage_of_contam_species = sorted_species[x][0]
        depth_of_contam_species = sorted_species[x][1]
        contam_species = sorted_species[x][2]
        if depth_of_contam_species > depth_of_top_species:
            warnings.append("warning: the top species hit (%s) has the highest %% coverage of all afanc species classifications (%s) and a median depth of %s, but a contaminating species (%s) - although with lower coverage (%s %%) - has higher depth (%s)" %(top_species, pc_coverage_of_top_species, depth_of_top_species, contam_species, pc_coverage_of_contam_species, depth_of_contam_species))

    # OTHER THAN THE TOP HIT, WHAT NON-HUMAN SPECIES ARE ALSO PRESENT IN THE SAMPLE, ACCORDING TO KRAKEN?
    no_of_human_reads = 0
    other_species = {}
    for key in kraken['Species']:
        species = key['name']
        taxid = key['taxon_id']
        reads = key['reads']
        if taxid == 9606: no_of_human_reads += reads
        # species = species.replace("Mycobacteriodes", "Mycobacterium") # Kraken sometimes uses "Mycobacteriodes" (e.g. Mycobacteriodes abscessus) whereas afanc uses "Mycobacterium" for the same. We need to standardise this to prevent downstream errors
        ## ignore any Kraken hits to mycobacterial species - they may be spurious. We will use only the mycobacterial classifications made by afanc
        if match_taxonomy(species): continue
        if taxid == 9606: continue # ignore human because we have a dedicated human read removal process elsewhere in the workflow
        if species != top_species: other_species[species] = taxid

    # OTHER THAN THE TOP HIT, WHAT NON-HUMAN SPECIES ARE ALSO PRESENT IN THE SAMPLE, ACCORDING TO AFANC?
    for species in afanc[sample_id]['phylogenetics']['species']:
        species = species.replace("_", " ")
        ## afanc does not assign a taxon ID to each species, so we will need to look this up. The taxon ID is the basis on which species' genomes are downloaded - we cannot proceed without it.
        if ((species not in tax_ids) & (species != top_species)):
            warnings.append("warning: unable to find a taxon ID for '%s', which means we will not be able to locate its genome, and thereby remove it as a contaminant. Check the Kraken report to see how this species has been reported" %species)
        if species not in tax_ids: continue
        taxid = tax_ids[species]
        if taxid == 9606: continue # ignore human because we have a dedicated human read removal process elsewhere in the workflow
        if species != top_species: other_species[species] = taxid

    # IDENTIFY GENOMES FOR EACH MEMBER OF THE NON-REDUNDANT LIST OF NON-TOP-HIT & NON-HUMAN SPECIES PRESENT IN THE SAMPLE
    species = []
    for spec in other_species:
        species.append(spec)
    sorted_species = species.copy()
    sorted_species.sort(reverse=True)

    ignored_mixed_myco = {}
    contaminant_genera = {}
    # list for printing urls into test_urllist.txt
    out_urls = []

    for spec in sorted_species:
        taxid = other_species[spec]
        if taxid == 9606: continue # ignore human reads once again
        if (taxid not in urls):
            warnings.append("warning: unable to find the latest RefSeq genome for taxon ID %s, and thereby remove it as a contaminant (the Kraken report assigns this taxon ID to species '%s')" %(taxid, spec))
        if (taxid not in urls): continue

        # for each contaminating species, we will identify the set of genomes available in RefSeq
	    # if one or more 'reference genomes' for this species are available and they are 'complete' (thereby implicitly having 'full' representation), we will use these
	    # else, we will use all complete genome available
	    # else, we will use all genomes that are available (which may simply be contigs) and warn the user that while we can proceed with contaminant removal using them, we will have reduced confidence in detecting reads from this species
        arr = urls[taxid]

        reference_genomes = {}
        for x in range(0, len(arr)):
            refseq_species = arr[x][0]
            infraspecific_name = arr[x][1]
            ftp_path = arr[x][2]
            assembly_level = arr[x][3]
            genome_rep = arr[x][4]
            refseq_category = arr[x][5]
            if ((assembly_level == 'Complete Genome') & (genome_rep == 'Full') & (refseq_species == 'reference genome')):
                if ftp_path not in reference_genomes: reference_genomes[ftp_path] = 0
                reference_genomes[ftp_path] += 1

        complete_genomes = {}
        for x in range(0, len(arr)):
            refseq_species = arr[x][0]
            infraspecific_name = arr[x][1]
            ftp_path = arr[x][2]
            assembly_level = arr[x][3]
            genome_rep = arr[x][4]
            refseq_category = arr[x][5]
            # checkpoint: if we have found reference genomes (scalar keys %refernece_genomes != 0) then we will ignore any genomes that are not (i.e. not in the %reference_genomes hash)
            if ((len(reference_genomes) > 0) & (ftp_path not in reference_genomes)): continue
            if ((assembly_level == 'Complete Genome') & (genome_rep == 'Full')):
                if ftp_path not in complete_genomes: complete_genomes[ftp_path] = 0
                complete_genomes[ftp_path] += 1

        for x in range(0, len(arr)):
            refseq_species = arr[x][0]
            infraspecific_name = arr[x][1]
            ftp_path = arr[x][2]
            assembly_level = arr[x][3]
            genome_rep = arr[x][4]
            refseq_category = arr[x][5]

            filename = ''
            re_filename = re.findall(r"^.+\/(.*?)$", ftp_path)
            if len(re_filename) > 0: filename = re_filename[0]
            if filename == '':
                warnings.append("warning: unable to parse FTP path to the FASTA of a contaminant species genome ('%s'), this being necessary to download it: %s" %(spec, ftp_path))
                continue
            full_path = ftp_path + '/' + filename + '_genomic.fna.gz'
            # checkpoint: if we have found complete genomes (scalar keys %complete_genomes != 0) then we will ignore any genomes that are incomplete (i.e. not in the %complete_genomes hash)
            if ((len(complete_genomes) > 0) & (ftp_path not in complete_genomes)): continue
            hash = {
                'infraspecific_name': infraspecific_name,
                'url': full_path,
                'assembly_level': assembly_level,
                'genome_representation': genome_rep,
                'refseq_category': refseq_category,
                'taxid': taxid
            }
            contaminant_genus = ''
            contaminant_species = ''
            re_species = re.findall("^(.*?) (.+?)$", spec)
            if len(re_species[0]) > 1:
                 contaminant_genus = re_species[0][0]
                 contaminant_species = re_species[0][1]
            if ((unmix_myco == 'no') & (match_taxonomy(top_species)) & (match_taxonomy(spec))):
                if spec not in ignored_mixed_myco: ignored_mixed_myco[spec] = 0
                ignored_mixed_myco[spec] += 1
            else:
                if contaminant_genus not in contaminant_genera: contaminant_genera[contaminant_genus] = {}
                if contaminant_species not in contaminant_genera[contaminant_genus]: contaminant_genera[contaminant_genus][contaminant_species] = 0
                contaminant_genera[contaminant_genus][contaminant_species] += 1

                if 'contaminants' not in out: out['contaminants'] = {}
                if spec not in out['contaminants']: out['contaminants'][spec] = []
                out['contaminants'][spec].append(hash)

                out_urls.append(full_path)

        if (len(complete_genomes) == 0):
            warnings.append("warning: no complete genome was found for the contaminant species '%s'. Alignment-based read removal will not necessarily detect all reads from this species" %spec)

    # ARE WE DELIBERATELY IGNORING MIXED-MYCOBACTERIAL CONTENT, CONSIDERING CLASSIFICATIONS TO OTHER MYCOBACTERIAL SPECIES TO *NOT* BE CONTAMINANTS?
    for name in ignored_mixed_myco:
        warnings.append("warning: sample contains a mixture of mycobacteria, its principal species (%s) mixed with reads from %s. As you have chosen to ignore this (--unmix_myco no), %s is NOT being counted as a contaminant" %(top_species, name, name))

    # HOW MANY CONTAMINATING GENERA ARE THERE?
    no_of_contaminant_genera = len(contaminant_genera)

    # DO WE DETECT ONE OR MORE CONTAMINANTS?
    if 'summary_questions' not in out: out['summary_questions'] = {}
    if 'contaminants' in out:
        out['summary_questions']['are_there_contaminants'] = 'yes'
    else:
        out['summary_questions']['are_there_contaminants'] = 'no'

    # *DID* WE DETECT ONE OR MORE CONTAMINANTS, WHEN (IF) WE RAN THIS SCRIPT BEFORE? NOTE THAT THIS SCRIPT CAN BE CALLED UP TO TWICE IN THE WORKFLOW.
    if prev_species != '':
        contam_found = 0
        if 'contaminants' in prev_species:
            for contam_species in prev_species['contaminants']:
                contam_found += 1
                if 'contaminants_previously_removed' not in out['summary_questions']: out['summary_questions']['contaminants_previously_removed'] = []
                out['summary_questions']['contaminants_previously_removed'].append(contam_species)
        if contam_found > 0:
            out['summary_questions']['were_contaminants_removed'] = 'yes'
        elif contam_found == 0:
            out['summary_questions']['were_contaminants_removed'] = 'no'
    else:
        out['summary_questions']['were_contaminants_removed'] = 'no'

    # IS THE TOP SPECIES HIT ONE OF THE 10 ACCEPTABLE POSSIBILITIES? IF SO, PROVIDE A LINK TO THE REFERENCE GENOME
    re_top_species = re.findall(r"^(Mycobact|Mycolicibac)\w+ (abscessus|africanum|avium|bovis|chelonae|chimaera|fortuitum|intracellulare|kansasii|tuberculosis).*?$", top_species)
    if len(re_top_species) > 0:
        identified_species = re_top_species[0][1]
        if supposed_species == 'null':
            out['summary_questions']['is_the_top_species_appropriate'] = 'yes'
        elif ((supposed_species != 'null') & (supposed_species == identified_species)):
            out['summary_questions']['is_the_top_species_appropriate'] = 'yes'
        elif ((supposed_species != 'null') & (supposed_species != identified_species)):
            warnings.append("warning: the top species hit is %s, contrary to the expectation: %s" %(identified_species, supposed_species))
            out['summary_questions']['is_the_top_species_appropriate'] = 'no'

        if out['summary_questions']['is_the_top_species_appropriate'] == 'yes':
            if out['summary_questions']['are_there_contaminants'] == 'yes':
                if no_of_contaminant_genera == 1:
                    contaminating_genus = ''
                    for genus in contaminant_genera: # !??
                        contaminating_genus = genus
                    if match_taxonomy(contaminant_genus):
                        contam_species = []
                        for contam_spec in contaminant_genera[contaminating_genus]:
                            contam_species.append(contam_spec)
                        sorted_contam_species = contam_species.copy()
                        sorted_contam_species.sort(reverse=True)
                        contam_species = ", ".join(sorted_contam_species)
                        num_contam_species = len(sorted_contam_species)
                        warnings.append("warning: sample contains a mixture of mycobacteria, its principal species (%s) mixed with reads from %d other Mycobacterium: %s" %(identified_species, num_contam_species, contam_species))
                    else:
                        warnings.append("warning: the top species hit (%s) is supported, but the sample shows signs of contamination with species from one other genus (%s)" %(identified_species, contaminating_genus))
                elif no_of_contaminant_genera > 1:
                    contam_genus = []
                    contam_with_myco = 0
                    for genus in contaminant_genera:
                        if match_taxonomy(genus):
                            contam_with_myco += 1
                            continue
                        contam_genus.append(genus)
                    sorted_contam_genus = contam_genus.copy()
                    sorted_contam_genus.sort(reverse=True)
                    non_myco_contam_genus = ", ".join(sorted_contam_genus)
                    num_non_myco_contam_genus = len(sorted_contam_genus)
                    if contam_with_myco == 0:
                        if num_non_myco_contam_genus == 1:
                            warnings.append("warning: the top species hit (%s) is supported, but the sample shows signs of contamination with species from %d other (non-mycobacterial) genus: %s" %(identified_species, num_non_myco_contam_genus, non_myco_contam_genus))
                        elif num_non_myco_contam_genus > 1:
                            warnings.append("warning: the top species hit (%s) is supported, but the sample shows signs of contamination with species from %d other (non-mycobacterial) genera: %s" %(identified_species, num_non_myco_contam_genus, non_myco_contam_genus))
                    elif contam_with_myco > 0:
                        if num_non_myco_contam_genus == 1:
                            warnings.append("warning: the top species hit (%s) is supported, but the sample contains multiple mycobacteria and shows signs of contamination with species from %d other (non-mycobacterial) genus: %s" %(identified_species, num_non_myco_contam_genus, non_myco_contam_genus))
                        elif num_non_myco_contam_genus > 1:
                            warnings.append("warning: the top species hit (%s) is supported, but the sample contains multiple mycobacteria and shows signs of contamination with species from %d other (non-mycobacterial) genera: %s" %(identified_species, num_non_myco_contam_genus, non_myco_contam_genus))
                else:
                    sys.exit("ERROR: sample is considered contaminated but the number of contaminant genera is %d" %no_of_contaminant_genera)

            ref_fa = os.path.join(myco_dir_path, identified_species + ".fasta")
            ref_dir = os.path.join(myco_dir_path, identified_species)
            ref_mmi = os.path.join(myco_dir_path, identified_species + ".mmi")

            if 'file_paths' not in out['top_hit']: out['top_hit']['file_paths'] = {}
            out['top_hit']['file_paths']['ref_fa'] = ref_fa
            out['top_hit']['file_paths']['clockwork_ref_dir'] = ref_dir

            ## commented out due to Gnomon presumably taking on drug resistance/susceptibility responsibility
            # susceptibility = afanc[sample_id]['susceptibility']
            # out['top_hit']['susceptibility'] = susceptibility

            out['top_hit']['phylogenetics'] = afanc[sample_id]['phylogenetics']

    else:
        out['summary_questions']['is_the_top_species_appropriate'] = 'no'

    # IN THIS WORKFLOW, AFANC WOULD ONLY BE CALLED IF KRAKEN CLASSIFIED > 100k READS AS MYCOBACTERIACEAE, SO FOR AFANC TO MAKE *NO CLASSIFICATION* SOMETHING SUSPECT IS GOING ON.
    # WHAT IS LIKELY TO HAVE HAPPENED IS THAT THE ALIGNMENT-BASED DECONTAMINATION PROCESS HAS TRIED TO DISAMBIGUATE A MIXTURE OF VERY SIMILAR MYCOBACTERIA AND INADVERTENTLY REMOVED TOO MANY READS. THERE WILL BE NOTHING SUBSTANTIVE LEFT FOR AFANC TO CLASSIFY.
    if ((num_afanc_species == afanc_finds_nothing) & (num_afanc_species == 1)):
        if out['summary_questions']['were_contaminants_removed'] == 'yes':
            warnings.append("warning: regardless of what Kraken reports, afanc did not make a species-level mycobacterial classification. If this is a mixed-mycobacterial sample, then an alignment-based contaminant-removal process may not be appropriate. Suggestion: re-run with --unmix_myco 'no'")
        elif out['summary_questions']['were_contaminants_removed'] == 'no':
            warnings.append("warning: regardless of what Kraken reports, afanc did not make a species-level mycobacterial classification")

    # IF THE TOP HIT IS APPROPRIATE AND THERE ARE NO CONTAMINANTS, WE CAN CONTINUE TO RUN CLOCKWORK
    if ((out['summary_questions']['is_the_top_species_appropriate'] == 'yes') & (out['summary_questions']['are_there_contaminants'] == 'no')):
        out['summary_questions']['continue_to_clockwork'] = 'yes'
    else:
        out['summary_questions']['continue_to_clockwork'] = 'no'

    # IF THE SAMPLE IS *NOT CURRENTLY CONSIDERED* CONTAMINATED BUT PREVIOUSLY WAS, WE SAY SO
    if ((out['summary_questions']['are_there_contaminants'] == 'no') & (out['summary_questions']['were_contaminants_removed'] == 'yes')):
        warnings.append("warning: the sample showed signs of contamination, although the contaminating reads were considered successfully removed")

    # IF NO WARNINGS HAVE BEEN RAISED, WE SAY SO
    if len(warnings) == 0:
        warnings.append("no warnings raised")
    out['warnings'] = warnings

    return out, out_urls

# call main function
if __name__ == "__main__":
    # set command line arguments
    description = 'This script will parse the Kraken and afanc output JSONs to identify (a) the dominant mycobacterial species in the sample (based on highest %% coverage, determined by afanc), and (b) all other species, which are considered contaminants\n'
    description += 'The output will be one JSON and one txt file, both produced in rundir. The JSON will state the afanc-determined dominant species in the sample, as well as list each contaminant species, showing the taxon IDs, fasta URL, and both the \'assembly level\' and \'genome representation\' for each\n'
    description += "The text file just contains the URLs for each contaminant fasta\n"
    description += "For each species considered a contaminant, the latest RefSeq genomes are obtained as follows:\n"
    description += "If available, obtain all NCBI 'reference genomes' of that species, provided they are defined as 'complete'\n"
    description += "Else: all 'complete genomes' of that species, regardless of whether they are the reference\n"
    description += "Else: any genome of that species, and we warn that it may not be complete (which reduces confidence in contaminant removal)\n"
    description += "A 'reference genome' is a manually-selected community standard for that species. Note that some prokaryotes can have more than one reference genome\n"
    description += "[species] refers to what you believe this sample to be. You will be warned if this differs from the Kraken/afanc predictions\n"
    description += "By defining [species] you will automatically select this to be the genome against which reads will be aligned using Clockwork\n"
    description += "[unmix myco] is either 'yes' or 'no', given in response to the question: do you want to disambiguate mixed-mycobacterial samples by read alignment?\n"
    description += "If 'no', any contaminating mycobacteria will be recorded but NOT acted upon\n"
    usage = "python identify_tophit_and_contaminants2.py [path to afanc JSON] [path to Kraken JSON] [path to RefSeq assembly summary file] [species] [unmix myco] [directory containing mycobacterial reference genomes] [aws_config]\n"
    usage += "E.G.:\tpython identify_tophit_and_contaminants2.py afanc_report.json afanc_report.json assembly_summary_refseq.txt 1 tuberculosis yes myco_dir\n\n\n"

    parser = argparse.ArgumentParser(description=description, usage=usage, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('afanc_json', metavar='afanc_json', type=str, help='Path to afanc json report')
    parser.add_argument('kraken_json', metavar='kraken_json', type=str, help='Path to Kraken json report')
    parser.add_argument('assembly_file', metavar='assembly_file', type=str, help='Path to RefSeq assembly summary file')
    parser.add_argument('species', metavar='species', type=str, help='Refers to what you believe this sample to be. You will be warned if this differs from the Kraken/afanc predictions')
    parser.add_argument('unmix_myco', metavar='unmix_myco', type=str, help='Is either \'yes\' or \'no\', given in response to the question: do you want to disambiguate mixed-mycobacterial samples by read alignment?\nIf \'no\', any contaminating mycobacteria will be recorded but NOT acted upon')
    parser.add_argument('myco_dir', metavar='myco_dir', type=str, help='Path to myco directory')
    parser.add_argument('prev_species_json', metavar='prev_species_json', type=str, help='Path to previous species json file. Can be set to \'null\'')
    args = parser.parse_args()

    # REQUIREMENTS
    sample_id = process_requirements(sys.argv)
    afanc_json = sys.argv[1]
    kraken_json = sys.argv[2]
    assembly_file = sys.argv[3]
    supposed_species = sys.argv[4]
    unmix_myco = sys.argv[5]
    myco_dir = sys.argv[6]
    prev_species_json = sys.argv[7]

    # read assembly summary
    urls, tax_ids = read_assembly_summary(assembly_file)

    # process reports
    out, out_urls = process_reports(afanc_json, kraken_json, supposed_species, unmix_myco, myco_dir, prev_species_json, urls, tax_ids, sample_id)

    # print urls into {sample_id}_urllist.txt
    out_file1 = sample_id + '_urllist.txt'
    with open(out_file1, 'w') as f:
        for url in out_urls:
            f.write(url + "\n")

    # print final file
    out_file2 = sample_id + '_species_in_sample.json'
    with open(out_file2, 'w') as f:
        json.dump(out, f, indent = 4)
