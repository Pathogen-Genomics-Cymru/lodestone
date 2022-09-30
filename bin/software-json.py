#!/usr/bin/env python3

import re
import json
import os
import glob
import sys


def go(path):

    software = []
    software_dict={}
    all_software_dict={}
    stop =  ['export', '\n']

    for filename in glob.glob(os.path.join(path, "Singularity.*")):
        extension = filename.split('.', 1)[1]
        version = filename.split('-')[-1]
        with open(os.path.join(filename), 'r') as infile:
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


    software_vers = 'lodestone-' + version
    all_software_dict = {software_vers : all_software_dict}

    with open('version.json', 'w') as outfile:
        json.dump(all_software_dict, outfile, indent=2, separators=(',', ' : '))


if __name__ == "__main__":
    go(sys.argv[1])

