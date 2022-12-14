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

        for event in data["results"]["Detection_events"]:
            if "closest_variant" in event:
                event = event["closest_variant"]

            id = event["name"].replace(" ", "_")
            id = id.replace(".", "")
            species_dict[id] = { "percent_coverage" : event["reference_cov"]*100, "median_depth" : event["median_DOC"] }

    new_json["arguments"] = data["arguments"]

    with open(f"{sample_id}_afanc_report.json", "w") as fout:
        json.dump(new_json, fout, indent = 4)

def main(argv):
    afanc_report =  argv[1]
    reformat_json(afanc_report)

if __name__=="__main__":
    main(sys.argv)
