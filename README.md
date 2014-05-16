Cluster Flow
============

Cluster Flow is a pipelining tool to automate and standardise bioinformatics analyses on high-performance cluster environments.

It is designed to be easy to use, quick to set up and flexible to configure.

Benefits of using Cluster Flow:
* Routine analyses are very quick to run, for example: `cf --genome GRCh37 fastq_bowtie *fq.gz`
* Pipelines use identical parameters, standardising analysis and making results more reproducable
* All commands and output is logged in files for future reference
* Intuitive commands and a comprehensive manual make Cluster Flow easy to use
* Works out of the box (almost - see the [video tutorial](http://youtu.be/b2g_zQiz9ys))

How Cluster Flow differs from other pipeline tools:
* Very lightweight and flexible
* Pipelines and configurations can easily be generated on a project-specific basis if required
* New modules and pipelines are trivial to write (see the [video tutorial](http://youtu.be/aBHOcsA2M6w))

There is a comprehensive [manual](http://www.bioinformatics.babraham.ac.uk/projects/clusterflow/Cluster_Flow_Manual.pdf) available for Cluster Flow describing installation, usage and customisation. A short version is below:

How Cluster Flow Works
----------------------
Cluster Flow is comprised of several layers: 
* cf 
	* The main cluster flow command. This is called to initiate a new pipeline run. 
* Pipelines 
	* Protocols that describe a series of modules to be run, along with any parameters 
* Modules 
	* The instructions for an individual task. These can be written in any language but must conform to a common API, described within this document. 
* Runs 
	* Created from the pipeline template for each file. Specifies configuration variables and traces output filenames. 
Cluster Flow will set off multiple queued jobs on the cluster with queue dependencies as defined in the pipeline. 


Installation
------------
1. Download Cluster Flow
	* You can either clone this github repository or download the files
	* Stable releases can be found on the releases page
2. Unpack and set up with environment modules (optional)
	* Cluster Flow is designed to work with environment modules, if available.
	* If you don't use the environment module system, add the Cluster Flow directory to your `PATH` so that the cf executable can be found.
	* If using environment modules, run `module load clusterflow` to load the module
3. Configure Cluster Flow (global)
	* Rename `clusterflow.config.example` and `genomes.config.example` so that they don't end in `.example`
	* Edit as necessary. See comments within the files.
4. Confugure Cluster Flow (user)
	* Run `cf --make_config` to set up Cluster Flow with your e-mail address and preferences
5. Add some reference genomes
	* Run `cf --add_genome` to specify the locations of your reference genomes

Typical Usage
-------------
Example command:

	cf --genome NCBIM37 sra_bowtie *sra

This calls Cluster Flow (`cf`), specifies the mouse reference genome (`--genome NCBIM37`), tells Cluster Flow to use the SRA input - bowtie alignment pipeline (`sra_bowtie`) and specifies the input files (`*sra`)

To you can see all availble pipelines, modules and refernece genomes with the following commands:

	cf --list_pipelines
	cf --list_modules
	cf --list_genomes

Command Line Parameters
-----------------------
For a full description of what each command line option does, see the [manual](http://www.bioinformatics.babraham.ac.uk/projects/clusterflow/Cluster_Flow_Manual.pdf) 

Flag | Description
---- | -----------
`--genome <ID>` | ID of a genome referred to in clusterflow.config 
`--genome_path <path>` | Path to a genome to be used for alignment 
`--bowtie_path <path>` | Path to a bowtie index base to be used for alignment 
`--gtf_path <path>` | Path to a GTF file to be used for alignment (eg. for Tophat) 
`--paired` | Force paired-end mode 
`--single` | Force single-end mode 
`--no_fn_check` | Disable input file type checking 
`--file_list` | Text file containing input files or download URLs 
`--params` | Specify extra module parameters for this run 
`--split_files <num>` | Create one run per `<num>` files 
`--max_runs <num>` | Divide input files into `<num>` runs. Set as 0 to disable. 
`--email <email>` | Set the e-mail address for notifications 
`--priority <num>` | Set the queue priority for cluster jobs 
`--cores <num>` | Set the maximum number of cores to use for all runs 
`--mem <string>` | Set the maximum memory to use for all runs 
`--notifications [cresa]` | Specify desired notifications 
`--list_pipelines` | Print available pipelines 
`--list_modules` | Print available modules 
`--list_genomes` | Print available genomes 
`--dry_run` | Prints jobs to terminal instead of submitting them to the cluster 
`--qstat` | Displays formatted qstat output of your jobs 
`--qstatall` | Displays formatted qstat output of all jobs  
`--qstatcols` | Colours output from --qstat or --qstatall 
`--qdel <id>` | Delete all jobs from a running pipeline. `<id>` is printed with `--qstat`
`--make_config` | Interactive prompt to generate a personalised CF config file 
`--add_genome` | Interactive wizard to add new genomes to your genomes.config files 
`--version` | Print version of Cluster Flow installed 
`--check_updates` | Look for available Cluster Flow updates 
`--help` | Print help message 

Credits
-------
Cluster Flow was written by [Phil Ewels](http://phil.ewels.co.uk) whilst working in the [Babraham Bioinformatics](http://www.bioinformatics.babraham.ac.uk/) group in Cambridge, UK. He now maintains it whilst working at the [Science for Life Laboratory](http://www.scilifelab.se/) in Stockholm, Sweden.
