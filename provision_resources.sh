#!/bin/bash
set -x

S3_URI=$1

# remove trailing /
S3_URI=`echo $S3_URI | sed 's:/$::'`

KRAKEN_DB='k2_pluspf_16gb_20220607.tar.gz'
KRAKEN_DB_URL='https://genome-idx.s3.amazonaws.com/kraken'
KRAKEN_DB_S3_URL="$S3_URI/kraken_db/k2_pluspf_16gb_20220607/"

BOWTIE2_INDEX='hg19_1kgmaj_snvs_bt2.zip'
BOWTIE2_INDEX_URL='https://genome-idx.s3.amazonaws.com/bt'
BOWTIE2_INDEX_S3_URI="$S3_URI/bowtie2_index/hg19_1kgmaj/"

TB_PIPELINE_RESOURCES_S3_URI="$S3_URI/tbpipeline/"

TMPDIR=$(mktemp -d)

# download and extract KrakenDB files
wget -c $KRAKEN_DB_URL/$KRAKEN_DB
mkdir $TMPDIR/k2db
tar xzf $KRAKEN_DB -C $TMPDIR/k2db

# upload KrakenDB files to S3 and remove local copy
aws s3 sync $TMPDIR/k2db $KRAKEN_DB_S3_URL
rm -r $TMPDIR/k2db
rm $KRAKEN_DB

# download and extract Bowtie index files
wget -c $BOWTIE2_INDEX_URL/$BOWTIE2_INDEX
mkdir $TMPDIR/bowtiedb
unzip $BOWTIE2_INDEX -d $TMPDIR/bowtiedb

# upload Bowtie index files to S3 and remove local copy
aws s3 sync $TMPDIR/bowtiedb $BOWTIE2_INDEX_S3_URI
rm -r $TMPDIR/bowtiedb
rm $BOWTIE2_INDEX

# move resources out of project directory and uplaod to S3
aws s3 sync resources $TB_PIPELINE_RESOURCES_S3_URI
mv resources ../tb-pipeline-resources
