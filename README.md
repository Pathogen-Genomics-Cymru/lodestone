# TB Pipeline #
  
Cleans and QCs reads with fastp and FastQC, classifies with Kraken2 & Mykrobe, removes non-bacterial content, and - by alignment to any minority genomes - disambiguates mixtures of bacterial reads.

Takes as input one directory containing pairs of fastq(.gz) or bam files.
Produces as output one directory per sample, containing the relevant reports & a pair of cleaned fastqs.

## Quick Start ## 
Requires NXF_VER>=20.11.0-edge

The workflow is designed to run with either docker `-profile docker` or singularity `-profile singularity`. Before running the workflow, the images will need to be built by running either `docker/docker_build.sh` or  `singularity/singularity_build.sh` 

E.g. to run the workflow:
```
nextflow run main.nf -profile singularity --filetype fastq --input_dir fq_dir --pattern "*_R{1,2}.fastq.gz" --unmix_myco yes \
--output_dir . --kraken_db /path/to/database --bowtie2_index /path/to/index --bowtie_index_name hg19_1kgmaj

nextflow run main.nf -profile docker --filetype bam --input_dir bam_dir --unmix_myco no \
--output_dir . --kraken_db /path/to/database --bowtie2_index /path/to/index --bowtie_index_name hg19_1kgmaj
```

## Params ##
The following parameters should be set in `nextflow.config` or specified on the command line:

* **input_dir**<br /> 
Directory containing fastq OR bam files
* **filetype**<br />
File type in input_dir. Either "fastq" or "bam"
* **pattern**<br />
Regex to match fastq files in input_dir, e.g. "*_R{1,2}.fq.gz"
* **output_dir**<br />
Output directory
* **unmix_myco**<br />
Do you want to disambiguate mixed-mycobacterial samples by read alignment? Either "yes" or "no"
* **species**<br />
Principal species in each sample, assuming genus Mycobacterium. Default 'null'. If parameter used, takes 1 of 10 values: abscessus, africanum, avium, bovis, chelonae, chimaera, fortuitum, intracellulare, kansasii, tuberculosis
* **kraken_db**<br />
Directory containing `*.k2d` Kraken2 database files (k2_pluspf_16gb_20200919 recommended, obtain from https://benlangmead.github.io/aws-indexes/k2)
* **bowtie2_index**<br />
Directory containing Bowtie2 index (obtain from ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19_1kgmaj_bt2.zip). The specified path should NOT include the index name
* **bowtie_index_name**<br />
Name of the bowtie index, e.g. hg19_1kgmaj<br />
<br />

For more information on the parameters run `nextflow run main.nf --help`

## Stub-run ##
To test the stub run:
```
nextflow run main.nf -stub -config testing.config
```

## Checkpoints ##
Checkpoints used throughout this workflow to fail a sample/issue warnings:
 
 processes preprocessing_checkFqValidity or preprocessing_checkBamValidity
 1. (Fail) If sample does not pass fqtools 'validate' or samtools 'quickcheck', as appropriate.
 
 process preprocessing_countReads\
 2. (Fail) If sample contains < 100k pairs of raw reads.
 
 process preprocessing_fastp\
 3. (Fail) If sample contains < 100k pairs of cleaned reads, required to all be > 50bp (cleaning using fastp with --length_required 50 --average_qual 10 --low_complexity_filter --correction --cut_right --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20).

 process preprocessing_kraken2\
 4. (Fail) If the top family hit is not Mycobacteriaceae\
 5. (Fail) If there are fewer than 100k reads classified as Mycobacteriaceae \
 6. (Warn) If the top family classification is mycobacterial, but this is not consistent with top genus and species classifications\
 7. (Warn) If the top family is Mycobacteriaceae but no G1 (species complex) classifications meet minimum thresholds of > 5000 reads or > 0.5% of the total reads (this is not necessarily a concern as not all mycobacteria have a taxonomic classification at this rank) \
 8. (Warn) If sample is mixed or contaminated - defined as containing reads > the 5000/0.5% thresholds from multiple non-human species\
 9. (Warn) If sample contains multiple classifications to mycobacterial species complexes, each meeting the > 5000/0.5% thresholds\
 10. (Warn) If no species classification meets the 5000/0.5% thresholds\
 11. (Warn) If no genus classification meets the 5000/0.5% thresholds\
 12. (Fail) If no family classification meets the 5000/0.5% thresholds (redundant given point 5)
 
 process preprocessing_identifyBacterialContaminants\
 13. (Fail) If the sample is not contaminated and the top species hit is not one of the 10 supported Mycobacteria:\ abscessus|africanum|avium|bovis|chelonae|chimaera|fortuitum|intracellulare|kansasii|tuberculosis\
 14. (Fail) If the sample is not contaminated and the top species hit is contrary to the species expected (e.g. "avium" rather than "tuberculosis" - only tested if you provide that expectation)\
 15. (Warn) If the top species hit is supported by < 75% coverage\
 16. (Warn) If the top species hit has a median coverage depth < 10-fold\
 17. (Warn) If we are unable to associate an NCBI taxon ID to any given contaminant species, which means we will not be able to locate its genome, and thereby remove it as a contaminant\
 18. (Warn) If we are unable to determine a URL for the latest RefSeq genome associated with a contaminant species' taxon ID\
 19. (Warn) If no complete genome could be found for a contaminant species. The workflow will proceed with alignment-based contaminant removal, but you're warned that there's reduced confidence in detecting reads from this species
 
 process preprocessing_downloadContamGenomes\
 20. (Fail) If a contaminant is detected but we are unable to download a representative genome, and thereby remove it
 
 process preprocessing_summarise\
 21. (Fail) If after having taken an alignment-based approach to decontamination, Kraken still detects a contaminant species\
 22. (Fail) If after having taken an alignment-based approach to decontamination, the top species hit is not one of the 10 supported Mycobacteria\
 23. (Fail) If, after successfully removing contaminants, the top species hit is contrary to the species expected (e.g. "avium" rather than "tuberculosis" - only tested if you provide that expectation)

