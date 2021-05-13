#!/bin/bash

container_list=("ppbedtools" "ppbowtie2" "ppbwa" "ppfastp" "ppfastqc" "ppfqtools" "ppkraken2" "ppmykrobe" "ppperljson")

for item in ${container_list[@]}; do
    str=${item}
    str2chr=${str:2}
    singularity pull ${str:0:2}${str2chr^}.sif docker://annacprice/${item}:latest 
done
