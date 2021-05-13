#!/bin/bash

container_list=("ppbedtools" "ppbowtie2" "ppbwa" "ppfastp" "ppfastqc" "ppfqtools" "ppkraken2" "ppmykrobe" "ppperljson")

for item in ${container_list[@]}; do
    docker pull annacprice/${item}:latest
done
