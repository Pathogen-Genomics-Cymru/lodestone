#!/bin/bash

KRAKEN_DB='k2_pluspf_16gb_20220607.tar.gz'
KRAKEN_DB_URL='https://genome-idx.s3.amazonaws.com/kraken'
KRAKEN_DB_S3_URI='s3://tbpipeline/kraken_db'

BOWTIE2_INDEX='hg19_1kgmaj'
BOWTIE2_INDEX_URL="https://genome-idx.s3.amazonaws.com/bt/${BOWTIE2_INDEX}_snvs_bt2.zip"
BOWTIE2_INDEX_S3_URI='s3://tbpipeline/bowtie2_index'

RESOURCE_S3_URI='s3://tbpipeline/resources/'

wget "${KRAKEN_DB_URL}/${KRAKEN_DB}"
tar xzf $KRAKEN_DB -C ./k2db
aws s3 cp ./k2db/* "$S3_URI/${KRAKEN_DB}"
rm -r ./k2db

wget "${BOWTIE2_INDEX_URL}/${BOWTIE2_INDEX}"
unzip ${BOWTIE2_INDEX}_snvs_bt2
aws s3 cp ./${BOWTIE2_INDEX}_snvs_bt2/* "$BOWTIE2_INDEX_S3_URI/${BOWTIE2_INDEX}"
rm -r ./${BOWTIE2_INDEX}_snvs_bt2

aws s3 cp ./resources/* $RESOURCE_S3_URI
