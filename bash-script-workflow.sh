#!/bin/bash

echo "Generating scenario inputs and manifests"
python ./tests/scrip-test/mainfest-script/manifest-automation.py
mv inputs.json inputs.json.template
mv MANIFEST.json MANIFEST.json.template

for i in {1..10}; do
   echo "Copying scenario $i" 
   cp ./tests/scrip-test/mainfest-script/output/inputs.scenario${i}.json ./inputs.json
   cp ./tests/scrip-test/mainfest-script/output/MANIFEST${i}.json ./MANIFEST.json

   echo "Submitting workflow $i" 
   #submit workflow here by calling the genetics cli
   agc workflow run tbpipeline --context ondemand
   echo "input.jsonscenario$i"
 done
FILE=inputs.json
if test -e "$FILE"; then
    echo "rm inputs.json"
fi

FILE2=MANIFEST.json
if test -e "$FILE2"; then 
    echo "rm MANIFEST.json"
fi
 
 # changes the manifest and input files from template to .json after all the test runs
 mv inputs.json.template inputs.json
 mv MANIFEST.json.template MANIFEST.json
 