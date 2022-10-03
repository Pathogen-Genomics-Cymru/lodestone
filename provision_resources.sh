#!/bin/bash

KRAKEN_DB='k2_pluspf_16gb_20220607.tar.gz'
KRAKEN_DB_URL='https://genome-idx.s3.amazonaws.com/kraken'
KRAKEN_DB_S3_URI='s3://tbpipeline/kraken_db'

BOWTIE2_INDEX=''
BOWTIE2_INDEX_URL=''
BOWTIE2_INDEX_S3_URI=''

wget "${KRAKEN_DB_URL}/${KRAKEN_DB}"
tar xzf $KRAKEN_DB -C ./k2db
aws s3 cp ./k3db/* "$S3_URI/${KRAKEN_DB}"

