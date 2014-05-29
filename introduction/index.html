---
title: Introduction
layout: toc
---

# What is Cluster Flow?
Cluster Flow is simple package to run pipelines in a cluster environment. It is comprised of several layers:

* `cf`
    * The main cluster flow command. This is called to initiate a new pipeline run.
* Pipelines
    * Protocols that describe a series of modules to be run, along with any parameters
* Modules
    * The instructions for an individual task. These can be written in any language but must conform to a common API, described within this document.
* Runs
    * Created from the pipeline template for each file. Specifies configuration variables and traces output filenames.

Cluster Flow will set off multiple queued jobs on the cluster with queue dependencies as defined in the pipeline.

# How does Cluster Flow work?
A typical Cluster Flow run will work as follows:

* Pipelines and modules are written for Cluster Flow
* The `cf` command is run, initiating a pipeline with a set of files
* CF decides on a number of runs, each with a set of input files or URLs
* Input filenames are checked to make sure that they all have the same file extension
* If input files are fastq, CF tries to work out whether they’re paired end or single end from the filenames. It dies with an error if it finds a mixture
* CF checks that the files exist (unless they’re URLs) and parses the pipeline file
* A run file is created for each run. This contains the configuration variables, a copy of the pipeline used and the starting filenames with the associated id `start_000`
* CF submits all of the jobs to the cluster
    * Each module in each run has its own cluster job
    * Jobs are submitted with dependencies so that they execute in the order specified by the pipeline
* Main CF program finishes and exits
* Modules execute, using the run file and the previous job ID to find their input files
* STDERR is appended to the log files.
    * Commands should print their commands with `###CFCMD` prepended
    * Important messages should start with `###CF`
* Upon the completion of each module, the filenames of any resulting output files are appended to the run file, along with the job ID so that the next module can find them
* Once all pipeline modules have finished, core cluster flow notification modules run
    * These clean up the output file and send e-mails

# Working Example
## Find what's available
Often the first step when running a Cluster Flow pipeline is to remind yourself what is possible. Cluster Flow has a handful of functions to help you do this.

### Available Pipelines
Running

	cf --list_pipelines

Will give a list of the pipelines installed in your copy of Cluster Flow, for example:

	================================
	Cluster Flow - available pipelines
	================================
	Installed pipelines:
	    Directory /clusterflow/0.1/pipelines/
		- bismark
		- bismark_pbat
		- bismark_singlecell
		- fastq_bismark
		- fastq_bowtie
		... (trimmed) ...

### Available Modules
Cluster Flow can run single modules as well as strung-together in pipelines. The analgous command is:

	cf --list_modules

Will give:

	================================
	Cluster Flow - available modules
	================================
	Available modules:
	    Directory /clusterflow/0.1/modules/
		- bismark_align
		- bismark_deduplicate
		- bismark_messy
		- bismark_methXtract
		- bismark_tidy
		- bowtie
		... (trimmed) ...

### Available Genomes
Most modules and pipelines require genomes. To see what is set up, run

	cf --list_genomes

Will give: 

```
================================
Cluster Flow - available genomes
================================

--------------------------------------------------
 /clusterflow/0.1/genomes.config
--------------------------------------------------

== Genome Paths ==
 Key                 Species             Assembly            Path
----------------------------------------------------------------------------------------------------
 GRCh37              Human               GRCh37         /Genomes/Human/GRCh37/
 GRCm38              Mouse               GRCm38         /Genomes/Mouse/GRCm38/
 NCBIM37             Mouse               NCBIM37        /Genomes/Mouse/NCBIM37/
 
== Bowtie Index Base Paths ==
 Key                 Species             Assembly            Path
----------------------------------------------------------------------------------------------------
 GRCh37              Human               GRCh37         /Genomes/Human/GRCh37/Homo_sapiens.GRCh37
 GRCm38              Mouse               GRCm38         /Genomes/Mouse/GRCm38/Mus_musculus.GRCm38
 NCBIM37             Mouse               NCBIM37        /Genomes/Mouse/NCBIM37/Mus_musculus.NCBIM37
 
== GTF File Paths ==
 Key                 Species             Assembly            Path
----------------------------------------------------------------------------------------------------
 GRCh37              Human               GRCh37         /Genomes/Human/GRCh37/Homo_sapiens.GRCh37.61.gtf
 GRCm38              Mouse               GRCm38         /Genomes/Mouse/GRCm38/Mus_musculus.GRCm38.70.gtf
 NCBIM37             Mouse               NCBIM37        /Genomes/Mouse/NCBIM37/Mus_musculus.NCBIM37.cleaned.gtf
```

## Running Cluster Flow
Once you're happy with what's available and what you want to run, you use the `cf` command.

For example, to run the `fastq_bowtie` pipeline (runs FastQC, trims the reads with Trim Galore! and then aligns with either bowtie 1 or bowtie 2, depending on the read length) on a selection of FastQ files, I would run:

	cf --genome GRCh37 fastq_bowtie *.fq.gz

A typical output might be:

	Processing 3 files in 12 runs. Submitting 1 file per run.

	Pipeline to be used:

	#fastqc
	#trim_galore
		#bowtie

	Processing file 1 - filename_1.fq.gz
	Processing file 2 - filename_2.fq.gz
	Processing file 3 - filename_3.fq.gz
	Your job 191764 ("cf_sra_bowtie_1401371996_sra_fqdump_579") has been submitted
	Your job 191765 ("cf_sra_bowtie_1401371996_fastqc_975") has been submitted
	Your job 191766 ("cf_sra_bowtie_1401371996_trim_galore_490") has been submitted
	Your job 191767 ("cf_sra_bowtie_1401371996_bowtie_576") has been submitted
	Your job 191768 ("cf_sra_bowtie_1401371996_email_run_complete_840") has been submitted
	Your job 191769 ("cf_sra_bowtie_1401371996_sra_fqdump_816") has been submitted
	Your job 191770 ("cf_sra_bowtie_1401371996_fastqc_484") has been submitted
	Your job 191771 ("cf_sra_bowtie_1401371996_trim_galore_930") has been submitted
	Your job 191772 ("cf_sra_bowtie_1401371996_bowtie_281") has been submitted
	Your job 191773 ("cf_sra_bowtie_1401371996_email_run_complete_911") has been submitted
	Your job 191774 ("cf_sra_bowtie_1401371996_sra_fqdump_940") has been submitted
	Your job 191775 ("cf_sra_bowtie_1401371996_fastqc_805") has been submitted
	Your job 191776 ("cf_sra_bowtie_1401371996_trim_galore_634") has been submitted
	Your job 191777 ("cf_sra_bowtie_1401371996_bowtie_109") has been submitted
	Your job 191778 ("cf_sra_bowtie_1401371996_email_run_complete_922") has been submitted
	Your job 191779 ("cf_sra_bowtie_1401371996_email_pipeline_complete_127") has been submitted

## Checking the status of a run
I can then check on the progress of my running job using:

	cf --qstat

This will give the following report (much easier to read than the default `qstat`)
	
	======================================================================
	 Cluster Flow Pipeline: sra_bowtie                                  
	 Submitted:             30 seconds ago                              
	 Working Directory:     /bi/group/bioinf/clusterflow/dev/test       
	 ID:                    sra_bowtie_1401371972                       
	======================================================================
	
	 -  sra_fqdump                         [1 core] running for 16s
	      - fastqc 
	      - trim_galore 
	           - bowtie 
	                - email_run_complete 
	
	 -  sra_fqdump                         [1 core]  [priority -500] queued for 27s
	      - trim_galore 
	           - bowtie 
	                - email_run_complete 
	                     - email_pipeline_complete 
	      - fastqc 
	
	 -  sra_fqdump                         [1 core]  [priority -500] queued for 27s
	      - fastqc 
	      - trim_galore 
	           - bowtie 
	                - email_run_complete
	

## Cancelling a pipeline
I might later realise that I messed something up with the input files, but since then I have started loads of other jobs. I don't want to pull out all of the job IDs myself manually, so instead I can use Cluster Flow to delete the jobs belonging to the ID printed with the command above:

	cf --qdel sra_bowtie_1401371972

Which might give the following outout:

	Deleting jobs from Pipeline id sra_bowtie_1401371972

	16 jobs deleted:
	ewelsp has registered the job 191748 for deletion
	ewelsp has deleted job 191753
	ewelsp has deleted job 191758
	ewelsp has deleted job 191749
	ewelsp has deleted job 191750
	ewelsp has deleted job 191751
	ewelsp has deleted job 191752
	ewelsp has deleted job 191754
	ewelsp has deleted job 191755
	ewelsp has deleted job 191756
	ewelsp has deleted job 191757
	ewelsp has deleted job 191759
	ewelsp has deleted job 191760
	ewelsp has deleted job 191761
	ewelsp has deleted job 191762
	ewelsp has deleted job 191763

