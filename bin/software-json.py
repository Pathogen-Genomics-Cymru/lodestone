#!/usr/bin/env python3

import re
import json
import os
import glob
import sys


def go(singpath, configpath):

    software = []
    database = []
    software_dict={}
    database_dict={}
    all_software_dict={}
    stop =  ['export', '\n']
    
    paths = list(set(glob.glob(os.path.join(singpath, "Singularity.*"))) - \
                    set(glob.glob(os.path.join(singpath, "*.img"))))

    for filename in paths:
        extension = filename.split('.', 1)[1]
        version = filename.split('-')[-1]
        with open(os.path.join(singpath, filename), 'r') as infile:
            copy = False
            for line in infile:
                if line.strip() == "%environment":
                    copy = True
                    continue
                elif line.strip() == "%runscript":
                    copy = False
                    continue
                elif copy:
                    software.append(line)

        r = re.compile('.*export.[A-Z]')
        for item in software[:]:
            if r.match(item):
                software.remove(item)

        software = [' '.join(word for word in phrase.split() if word not in stop) for phrase in software]
        software = [item.replace('=', ':') for item in software]

        software_dict = dict(item.split(':') for item in software)
        software_dict = {extension : software_dict}

        all_software_dict.update(software_dict)

        software.clear()
        software_dict.clear()

    database_list = ["afanc_myco_db", "kraken_db", "bowtie2_index", "bowtie_index_name", "amr_cat"]
    replaced = [' ', "'", '"', '\n']

    with open(configpath) as infile:
        for line in infile:
            for elem in replaced:
                line = line.replace(elem, '')
            for item in database_list:
                if line.startswith(item):
                    database.append(line)

    database = [item.replace('=', ':') for item in database]
    database_dict = dict(item.split(':', 1) for item in database)
    database_dict = {"databases" : database_dict}

    all_software_dict.update(database_dict)

    software_vers = 'lodestone-' + version
    all_software_dict = {software_vers : all_software_dict}

    with open('version.json', 'w') as outfile:
        json.dump(all_software_dict, outfile, indent=2, separators=(',', ' : '))


if __name__ == "__main__":
    go(sys.argv[1], sys.argv[2])

