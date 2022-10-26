#!/bin/bash
set -e

KRAKEN_DB='k2_pluspf_16gb_20220607.tar.gz'
KRAKEN_DB_URL='https://genome-idx.s3.amazonaws.com/kraken'
KRAKEN_DB_S3_URL='s3://tbpipeline/kraken_db/k2_pluspf_16gb_20220607/'

BOWTIE2_INDEX='hg19_1kgmaj_snvs_bt2.zip'
BOWTIE2_INDEX_URL='https://genome-idx.s3.amazonaws.com/bt'
BOWTIE2_INDEX_S3_URI='s3://tbpipeline/bowtie2_index/hg19_1kgmaj/'


wget -c "${KRAKEN_DB_URL}/${KRAKEN_DB}"
TMPDIR=$(mktemp -d) #runs the command and prints the name
mkdir $TMPDIR/k2db
tar xzf $KRAKEN_DB -C $TMPDIR/k2db
aws s3 sync $TMPDIR/k2db $KRAKEN_DB_S3_URL
rm -r $TMPDIR/k2db
rm $KRAKEN_DB

wget -c "${BOWTIE2_INDEX_URL}/${BOWTIE2_INDEX}"
mkdir $TMPDIR/bowtiedb
unzip $BOWTIE2_INDEX -d $TMPDIR/bowtiedb
aws s3 sync $TMPDIR/bowtiedb $BOWTIE2_INDEX_S3_URI
rm -r bowtiedb
rm $BOWTIE2_INDEX

aws s3 sync tb-pipeline-resources s3://tbpipeline/resources/
