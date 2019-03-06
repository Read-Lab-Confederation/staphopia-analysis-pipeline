[![DOI](https://zenodo.org/badge/54082104.svg)](https://zenodo.org/badge/latestdoi/54082104)

# Staphopia Analysis Pipeline
Staphopia's analysis pipeline is a set of open source tools wrapped in a Nextflow workflow. Once analysis is completed the results are stored into a database.

![Staphopia Workflow](/data/staphopia-workflow.png "Staphopia Workflow")

# Installation
Currently Staphopia is setup using a Docker container.

## Pull Docker Container
```
docker pull rpetit3/staphopia:20190306
```

## Local Genome Analysis with `staphopia.py`
Staphopia pipeline, `staphopia.py`, processes local FASTQ files.

### Usage
```
docker run rpetit3/staphopia:112017 staphopia.py --help
usage: staphopia.py [-h] --fq1 FASTQ --sample SAMPLE [--fq2 FASTQ]
                    [--coverage INT] [--cpu INT] [--is_miseq] [--resume]

A wrapper for executing Staphopia Nextflow workflow.

optional arguments:
  -h, --help       show this help message and exit
  --fq1 FASTQ      Input FASTQ file.
  --sample SAMPLE  Sample name of the input.
  --fq2 FASTQ      Second FASTQ file in paired end reads.
  --coverage INT   Coverage to subsample to (Default: 100x.)
  --cpu INT        Number of processors to use.
  --is_miseq       Input is Illumina MiSeq sequencing.
  --resume         Tell nextflow to resume the run.
```

### Example
```
rpetit@staphopia:~/JE2$ ls -lh
total 701M
-rw-rw-r-- 1 rpetit rpetit 360M Oct  5  2017 JE2_R1.fastq.gz
-rw-rw-r-- 1 rpetit rpetit 342M Oct  5  2017 JE2_R2.fastq.gz
rpetit@staphopia:~/JE2$ docker run --rm -v $PWD:/data rpetit3/staphopia:20190306 staphopia.py \
                                   --fq1 JE2_R1.fastq.gz --fq2 JE2_R2.fastq.gz --sample JE2 \
                                   --cpu 22 --is_miseq
rpetit@staphopia:~/JE2$ ls -lh
total 990M
drwxr-xr-x 6 rpetit rpetit 4.0K Feb 12 19:26 JE2
-rw-r--r-- 1 rpetit rpetit   45 Feb 12 19:27 JE2.md5
-rw-rw-r-- 1 rpetit rpetit 360M Feb 12 17:41 JE2_R1.fastq.gz
-rw-rw-r-- 1 rpetit rpetit 342M Feb 12 17:41 JE2_R2.fastq.gz
-rw-r--r-- 1 rpetit rpetit  38K Feb 12 19:27 JE2-staphopia.txt
-rw-r--r-- 1 rpetit rpetit 288M Feb 12 19:26 JE2.tar.gz
```

## ENA Genome Analysis with `staphopia-ena.py`
Staphopia includes a modified pipeline, `staphopia-ena.py`, that automates the download of FASTQ files from ENA. Users give an ENA Experiment accession as the input and any corresponding FASTQs are downloaded and processed. 

### Usage
```
docker run rpetit3/staphopia:112017 staphopia-ena.py --help
usage: staphopia-ena.py [-h] [--cpu INT] [--resume] EXPERIMENT_ACCESSION

A wrapper for executing Staphopia Nextflow workflow.

positional arguments:
  EXPERIMENT_ACCESSION  ENA experiment accession to process.

optional arguments:
  -h, --help            show this help message and exit
  --cpu INT             Number of processors to use.
  --resume              Tell nextflow to resume the run.
```
### Example
```
rpetit@staphopia:~/$ mkdir SRX1114352
rpetit@staphopia:~/SRX1114352$ cd SRX1114352
rpetit@staphopia:~/SRX1114352$ docker run --rm -v $PWD:/data rpetit3/staphopia:20190306 staphopia-ena.py SRX1114352 --cpu 22
rpetit@staphopia:~/SRX1114352$ ls -lh
total 204M
-rw-r--r-- 1 root root   52 Apr 16 17:46 SRX1114352.md5
-rw-r--r-- 1 root root  39K Apr 16 17:46 SRX1114352-staphopia.txt
-rw-r--r-- 1 root root 204M Apr 16 17:46 SRX1114352.tar.gz
```

# List of Tools
Below is a list of the software programs used in Staphopia's pipeline. We would like to thank the authors of these tools!

| Tool        | Version           |
| ------------- |-------------|
|[Ariba](https://github.com/sanger-pathogens/ariba)|2.10.2|
|[assembly-summary](https://github.com/rpetit3/assembly-summary)|0.1|
|[BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)|37.66|
|[bedtools](http://bedtools.readthedocs.org/en/latest/)|2.26|
|[BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)|2.7.1|
|[BWA](https://github.com/lh3/bwa)|0.7.17|
|[ena-dl](https://github.com/rpetit3/ena-dl)|0.1|
|[GATK](https://software.broadinstitute.org/gatk/)|3.8|
|[illumina-cleanup](https://github.com/rpetit3/illumina-cleanup)|0.3|
|[Jellyfish](http://www.genome.umd.edu/jellyfish.html)|2.2.6|
|[MentaLiST]()|0.1.3|
|[Nextflow](https://www.nextflow.io/)|0.28.2|
|[Picard](https://broadinstitute.github.io/picard/)|2.14.1|
|[PROKKA](https://github.com/tseemann/prokka/)|1.12|
|[Samtools](http://www.htslib.org/)|1.6|
|[SPAdes](http://bioinf.spbau.ru/spades)|3.11.1|
|[vcf-annotator](https://github.com/rpetit3/vcf-annotator)|0.4|

# Contact
Robert Petit robert.petit@emory.edu

