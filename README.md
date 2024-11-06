# Lodestone #
![Build Status](https://github.com/Pathogen-Genomics-Cymru/lodestone/workflows/build-push-quay/badge.svg)
![Build Status](https://github.com/Pathogen-Genomics-Cymru/lodestone/workflows/pytest/badge.svg)
![Build Status](https://github.com/Pathogen-Genomics-Cymru/lodestone/actions/workflows/stub-run.yml/badge.svg?branch=main)
## Table of Contents
- [What is Lodestone](#what-is-lodestone)
- [Quick Start](#quick-start)
- [Executors](#executors)
- [System Requirements](#system-requirements)
- [Parameters](#parameters)
- [Stub Runs](#stub-runs)
- [Checkpoints](#checkpoints)
- [Acknowledgments](#acknowledgements)
- [License](#-license)

## What is Lodestone?
 
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

## Paramaters ##
The following parameters should be set in `nextflow.config`. They can be accessed by `nextflow run main.nf --help`:

```
--input_dir                      [string]          Input directory containing FASTQs or BAMs
--pattern                        [string]          Glob pattern for FASTQs or BAM
--output_dir                     [string]          Output directory
--permissive                     [boolean]         Flag. If True, errors in decontamination will be demoted to warnings
--filetype                       [string]          Either "fastq" or "bam". Assumes FASTQs are PE Illumina reads and BAMs are mapped against one of the references in resources/  (accepted: bam, fastq) [default: fastq]
--unmix_myco                     [boolean]         Flag. If True then minority Mycobacteriaceae reads will be removed. If False, they will be discarded
--species                        [string]          Species which will be mapped against, corresponding to references in resources/: can be one of  abscessus, africanum, avium, bovis, chelonae, chimaera, fortuitum, intracellulare, kansasii, tuberculosis or null. If 'null' the top hit as determined by Afanc will be used  (accepted: null, abscessus, africanum,
avium, bovis, chelonae, chimaera, fortuitum, intracellulare, kansasii, tuberculosis)
--sing_dir                       [string]          Directory to singularity definition files. Used to parse versions for reporting [default: ${baseDir}/resources]
--config_file                    [string]          Path to Nextflow config file. Used for parsing arguments to write to results if needed [default: ${baseDir}/nextflow.config]
--help                           [boolean, string] Show the help message for all top level parameters. When a parameter is given to `--help`, the full help message of that parameter will be printed.
--helpFull                       [boolean]         Show the help message for all non-hidden parameters.
--showHidden                     [boolean]         Show all hidden parameters in the help message. This needs to be used in combination with `--help` or `--helpFull`.

resources
  --resource_dir                 [string] Path to resources directroy where utility files are stored [default: ${baseDir}/resources]
  --refseq                       [string] Path to NCBI refseq summary file [default: ${baseDir}/resources/assembly_summary_refseq.txt]

resistance
  --resistance_profiler          [string]  Tool used for tb-profiler. Either tb-profiler or tbtamr  (accepted: tb-profiler, tbtamr) [default: tb-profiler]
  --collate                      [boolean] Flag. If True resistance reports will be summarised

bowtie
  --bowtie_index                 [string] Bowtie index directory [default: ${baseDir}/bowtie2/]
  --bowtie_index_name            [string] Prefix for the Bowtie2 index (minus the file extensions). [default: hg19_1kgmaj]

afanc
  --afanc_percent_threshold      [number]  Minimum percentage threshold for reads in order for a taxa to be considered in Afanc if the pipeline has failed earlier on (for reporting) [default: 5]
  --afanc_n_reads_threshold      [integer] Minimum reads threshold for reads in order for a taxa to be considered in Afanc [default: 500]
  --afanc_fail_percent_threshold [number]  Minimum percentage threshold for reads in order for a taxa to be considered in Afanc [default: 2]
  --afanc_fail_n_reads_threshold [integer] Minimum reads threshold for reads in order for a taxa to be considered in Afanc if the pipeline has failed earlier on (for reporting) [default: 200]

kraken
  --kraken_percent_threshold     [number]  Percentage threshold of reads required for taxa to be included in Kraken reports [default: 10]
  --kraken_n_reads_threshold     [integer] Raw reads threshold required for taxa to be included in Kraken reports [default: 10000]
  --kraken_db                    [string]  Kraken2 database path [default: kraken2/]
```

## Stub runs ##
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

## License
The tool is licensed under the V3 GNU Affero GPL license. Please see the [LICENSE](LICENSE) file for more details.
