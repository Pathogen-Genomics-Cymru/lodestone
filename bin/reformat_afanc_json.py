#!/usr/bin/env python3

""" Reformats the Afanc report json for consumption by identify_tophit_and_contaminants
"""

from collections import defaultdict
import json
import sys

def reformat_json(afanc_report):

    new_json = defaultdict(dict)
    sample_id = afanc_report.split(".json")[0].split("/")[-1]

    with open(afanc_report, 'r') as fin:
        data = json.load(fin)

        new_json[sample_id]["phylogenetics"] = {}
        new_json[sample_id]["phylogenetics"]["species"] = {}
        species_dict = new_json[sample_id]["phylogenetics"]["species"]

        cr = data["results"]["Detection_events"]["Clustering_results"]
        vp = data["results"]["Detection_events"]["Variant_profile"]

        for c_event in cr:
            event = c_event
            species = "_".join(c_event["name"].split(" ")[:2])
            id  = species.replace(".", "")

            ## capture closest variant info
            if "closest_variant" in c_event:
                event = c_event["closest_variant"]
                id = event["name"].replace(" ", "_")

            ## check to see if the top species hit exists within the variant profile
            if species in vp and len(vp[species]) > 0:
                ## if so, construct a new hit id
                id  = species + "_" + sorted(list(vp[species].keys()))[-1]

            species_dict[id] = { "percent_coverage" : event["reference_cov"]*100, "median_depth" : event["median_DOC"] }

    new_json["arguments"] = data["arguments"]
    new_json["versions"] = data["versions"]

    with open(f"{sample_id}_afanc_report.json", "w") as fout:
        json.dump(new_json, fout, indent = 4)

def main(argv):
    afanc_report =  argv[1]
    reformat_json(afanc_report)

if __name__=="__main__":
    main(sys.argv)
