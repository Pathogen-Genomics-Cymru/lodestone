#!/usr/bin/env python3

import os
import sys
import urllib.request
import json
from pathlib import Path
from vcfmix import lineageScan

def go(vcf_file):
    # create a lineagescan object
    v = lineageScan()

    # assuming postfix of ".minos.vcf"
    sampleid = vcf_file.replace(".minos.vcf", "")
    print(sampleid)

    res = v.parse(vcffile=vcf_file, sample_id=sampleid)

    # print details of the regions scanned
    print(v.region_stats)

    # export details of the regions scanned
    v.region_stats.to_csv(f"{sampleid}_vcfmix-regions.csv")

    # compute F2 and F47 statistics (see publication)
    summary1 = v.f_statistics()
    print(f"v.f_statistics(): {summary1}")

    # save F2 and F47 statistics to json
    with open(f"{sampleid}_f-stats.json", 'w') as f:
        f.write(json.dumps(summary1, indent=4))

if __name__ == "__main__":
    go(sys.argv[1])

