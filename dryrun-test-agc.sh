#!/bin/bash

echo "Generating scenario inputs.json and MANIFEST.jsons"
python ./tests/scrip-test/mainfest-script/manifest-automation.py

# backup inputs.json and MANIFEST.json in project root directory
mv inputs.json inputs.json.template
mv MANIFEST.json MANIFEST.json.template

for i in {1..10}; do
    echo "Submitting workflow run for test scenario $i"

    # copy scenario input.json and MANIFEST.json
    cp ./tests/scrip-test/mainfest-script/output/inputs.scenario${i}.json ./inputs.json
    cp ./tests/scrip-test/mainfest-script/output/MANIFEST${i}.json ./MANIFEST.json

    # submit workflow run
    agc workflow run tbpipeline --context ondemand

done

# changes the manifest and input files from template to .json after all the test runs
mv inputs.json.template inputs.json
mv MANIFEST.json.template MANIFEST.json
