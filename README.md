# Lodestone #
![Build Status](https://github.com/Pathogen-Genomics-Cymru/lodestone/workflows/build-push-quay/badge.svg)
![Build Status](https://github.com/Pathogen-Genomics-Cymru/lodestone/workflows/pytest/badge.svg)
![Build Status](https://github.com/Pathogen-Genomics-Cymru/lodestone/workflows/stub-run/badge.svg)
  
This pipeline takes as input reads presumed to be from one of 10 mycobacterial genomes: abscessus, africanum, avium, bovis, chelonae, chimaera, fortuitum, intracellulare, kansasii, tuberculosis. Input should be in the form of one directory containing pairs of fastq(.gz) or bam files.

Pipeline cleans and QCs reads with fastp and FastQC, classifies with Kraken2 & Afanc, removes non-bacterial content, and - by alignment to any minority genomes - disambiguates mixtures of bacterial reads. Cleaned reads are aligned to either of the 10 supported genomes and variants called. Produces as output one directory per sample, containing cleaned fastqs, sorted, indexed BAM, VCF, F2 and F47 statistics, an antibiogram and summary reports.

Note that while Mykrobe is included within this pipeline, it runs as an independent process and is not used for any downstream reporting.

**WARNING**: There are currently known errors with vcfmix, as such `errorStrategy 'ignore'` has been added to the processes vcfpredict:vcfmix to stop the pipeline from crashing. Please check the stdout from nextflow to see whether these processes have ran successfully.

## Quick Start ## 
This is a Nextflow DSL2 pipeline, it requires a version of Nextflow that supports DSL2 and the stub-run feature. It is recommended to run the pipeline with  `NXF_VER=20.11.0-edge`, as the pipeline has been tested using this version. E.g. to download
```
export NXF_VER="20.11.0-edge"
curl -fsSL https://get.nextflow.io | bash
```

The workflow is designed to run with either docker `-profile docker` or singularity `-profile singularity`. The container images are pulled from quay.io and a singularity cache directory is set in the `nextflow.config`. 

*Please note, when running the pipeline with Singularity, either the ```$TMPDIR``` or ```$SINGULARITY_TMPDIR``` must be specified; e.g. ```export TMPDIR="./"``` to run in the working directory.*

E.g. to run the workflow:
```
NXF_VER=20.11.0-edge nextflow run main.nf -profile singularity --filetype fastq --input_dir fq_dir --pattern "*_R{1,2}.fastq.gz" --unmix_myco yes \
--output_dir . --kraken_db /path/to/database --bowtie2_index /path/to/index --bowtie_index_name hg19_1kgmaj

NXF_VER=20.11.0-edge nextflow run main.nf -profile docker --filetype bam --input_dir bam_dir --unmix_myco no \
--output_dir . --kraken_db /path/to/database --bowtie2_index /path/to/index --bowtie_index_name hg19_1kgmaj
```

There is also a pre-configured climb profile to run Lodestone on a CLIMB Jupyter Notebook Server. Add ```-profile climb``` to your command invocation. The input directory can point to an S3 bucket natively (e.g. ```--input_dir s3://my-team/bucket```). By default this will run the workflow in Docker containers and take advantage of kubernetes pods. The Kraken2, Bowtie2 and Afanc databases will by default point to the ```pluspf16```, ```hg19_1kgmaj_bt2``` and ```Mycobacteriaciae_DB_7.0``` directories by default. These are mounted on a public S3 bucket hosted on CLIMB.

### Executors ###

By default, the pipeline will just run on the local machine. To run on a cluster, modifications will have to be made to the `nextflow.config` to add in the executor. E.g. for a SLURM cluster add `process.executor = 'slurm'`. For more information on executor options see the Nextflow docs: https://www.nextflow.io/docs/latest/executor.html

### System Requirements ###
Minimum recommended requirements: 32GB RAM, 8CPU

## Params ##
The following parameters should be set in `nextflow.config` or specified on the command line:

* **input_dir**<br /> 
Directory containing fastq OR bam files
* **filetype**<br />
File type in input_dir. Either "fastq" or "bam"
* **pattern**<br />
Regex to match fastq files in input_dir, e.g. "*_R{1,2}.fq.gz". Only mandatory if --filetype is "fastq"
* **output_dir**<br />
Output directory for results
* **unmix_myco**<br />
Do you want to disambiguate mixed-mycobacterial samples by read alignment? Either "yes" or "no":
  * If "yes" workflow will remove reads mapping to any minority mycobacterial genomes but in doing so WILL ALMOST CERTAINLY ALSO reduce coverage of the principal species
  * If "no" then mixed-mycobacterial samples will be left alone. Mixtures of mycobacteria + non-mycobacteria will still be disambiguated
* **species**<br />
Principal species in each sample, assuming genus Mycobacterium. Default 'null'. If parameter used, takes 1 of 10 values: abscessus, africanum, avium, bovis, chelonae, chimaera, fortuitum, intracellulare, kansasii, tuberculosis. Using this parameter will apply an additional sanity test to your sample
  * If you DO NOT use this parameter (default option), pipeline will determine principal species from the reads and consider any other species a contaminant
  * If you DO use this parameter, pipeline will expect this to be the principal species. It will fail the sample if reads from this species are not actually the majority
* **kraken_db**<br />
Directory containing `*.k2d` Kraken2 database files (k2_pluspf_16gb recommended, obtain from https://benlangmead.github.io/aws-indexes/k2)
* **bowtie2_index**<br />
Directory containing Bowtie2 index (obtain from ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19_1kgmaj_bt2.zip). The specified path should NOT include the index name
* **bowtie_index_name**<br />
Name of the bowtie index, e.g. hg19_1kgmaj<br />
* **vcfmix**<br />
Run [vcfmix](https://github.com/AlexOrlek/VCFMIX), yes or no. Set to no for synthetic samples<br />
* **resistance_profiler**<br />
Run resistance profiling for Mycobacterium tubercuclosis. Either ["tb-profiler"](https://tbdr.lshtm.ac.uk/) or "none".
* **afanc_myco_db**<br />
Path to the [afanc](https://github.com/ArthurVM/Afanc) database used for speciation. Obtain from  https://s3.climb.ac.uk/microbial-bioin-sp3/Mycobacteriaciae_DB_7.0.tar.gz
* **update_tbprofiler**<br />
Update tb-profiler. Either "yes" or "no". "yes" may be useful when running outside of a container for the first time as we will not have constructed a tb-profiler database matching our reference. This is not needed with the climb, docker and singluarity profiles as the reference has already been added. Alternatively you can run ```tb-profiler update_tbdb --match_ref <lodestone_dir>/resources/tuberculosis.fasta```.
* **refseq**<br />
Path to assembly summary refseq file (taken from [here](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt)). A local version is stored for reproducibility purposes in ```resources/``` but for best results download the latest version.
* **refseq**<br />
One of "yes" or "no". If "yes", continue to clockwork flags will be ignored and alignment will be performed anyway. If there are not enough reads and/or not a reference found the programme will still exit.

For more information on the parameters run `nextflow run main.nf --help`

The path to the singularity images can also be changed in the singularity profile in `nextflow.config`. Default value is `${baseDir}/singularity`

## Stub-run ##
To test the stub run:
```
NXF_VER=20.11.0-edge nextflow run main.nf -stub -config testing.config
```

## Checkpoints ##
Checkpoints used throughout this workflow to fail a sample/issue warnings:

**processes preprocessing:checkFqValidity or preprocessing:checkBamValidity**
1. (*Fail*) If sample does not pass fqtools 'validate' or samtools 'quickcheck', as appropriate.
 
**process preprocessing:countReads**

2. (*Fail*) If sample contains < 100k pairs of raw reads.
 
**process preprocessing:fastp**

3. (*Fail*) If sample contains < 100k pairs of cleaned reads, required to all be > 50bp (cleaning using fastp with --length_required 50 --average_qual 10 --low_complexity_filter --correction --cut_right --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20).

**process preprocessing:kraken2**

4. (*Fail*) If the top family hit is not Mycobacteriaceae
5. (*Fail*) If there are fewer than 100k reads classified as Mycobacteriaceae 
6. (*Warn*) If the top family classification is mycobacterial, but this is not consistent with top genus and species classifications
7. (*Warn*) If the top family is Mycobacteriaceae but no G1 (species complex) classifications meet minimum thresholds of > 5000 reads or > 0.5% of the total reads (this is not necessarily a concern as not all mycobacteria have a taxonomic classification at this rank). In ```nextflow.config``` these files can be modified.
8. (*Warn*) If sample is mixed or contaminated - defined as containing reads > the 5000/0.5% thresholds from multiple non-human species
9. (*Warn*) If sample contains multiple classifications to mycobacterial species complexes, each meeting the > 5000/0.5% thresholds
10. (*Warn*) If no species classification meets the 5000/0.5% thresholds
11. (*Warn*) If no genus classification meets the 5000/0.5% thresholds
 
**process preprocessing:identifyBacterialContaminants**

12. (*Fail*) If regardless of what Kraken reports, Afanc does not make a species-level mycobacterial classification (note that we do not use Kraken mycobacterial classifications other than to determine whether 100k reads are family Mycobacteriaceae; for higher-resolution classification, we defer to Afanc)
13. (*Fail*) If the sample is not contaminated and the top species hit is not one of the 10 supported Mycobacteria: abscessus|africanum|avium|bovis|chelonae|chimaera|fortuitum|intracellulare|kansasii|tuberculosis
14. (*Fail*) If the sample is not contaminated and the top species hit is contrary to the species expected (e.g. "avium" rather than "tuberculosis" - only tested if you provide that expectation)
15. (*Warn*) If the top Afanc species hit, on the basis of highest % coverage, does not also have the highest median depth
16. (*Warn*) If we are unable to associate an NCBI taxon ID to any given contaminant species, which means we will not be able to locate its genome, and thereby remove it as a contaminant
17. (*Warn*) If we are unable to determine a URL for the latest RefSeq genome associated with a contaminant species' taxon ID
18. (*Warn*) If no complete genome could be found for a contaminant species. The workflow will proceed with alignment-based contaminant removal, but you're warned that there's reduced confidence in detecting reads from this species
 
**process preprocessing:downloadContamGenomes**

19. (*Fail*) If a contaminant is detected but we are unable to download a representative genome, and thereby remove it
 
**process preprocessing:summarise**

20. (*Fail*) If after having taken an alignment-based approach to decontamination, Kraken still detects a contaminant species
21. (*Fail*) If after having taken an alignment-based approach to decontamination, the top species hit is not one of the 10 supported Mycobacteria
22. (*Fail*) If, after successfully removing contaminants, the top species hit is contrary to the species expected (e.g. "avium" rather than "tuberculosis" - only tested if you provide that expectation)

**process clockwork:alignToRef**

23. (*Fail*) If < 100k reads could be aligned to the reference genome
24. (*Fail*) If, after aligning to the reference genome, the average read mapping quality < 10
25. (*Fail*) If < 50% of the reference genome was covered at 10-fold depth

**process clockwork:minos**

26. (*Warn*) If sample is not TB, then it is not passed to a resistance profiler

## Acknowledgements ##
For a list of direct authors of this pipeline, please see the contributors list. All of the software dependencies of this pipeline are recorded in the version.json

The preprocessing sub-workflow is based on the preprocessing nextflow DSL1 pipeline written by Stephen Bush, University of Oxford. The clockwork sub-workflow uses aspects of the variant calling workflow from https://github.com/iqbal-lab-org/clockwork, lead author Martin Hunt, Iqbal Lab at EMBL-EBI

