---
title: Installation Instructions
layout: toc
---

# Installation Overview
The heart of Cluster Flow is a simple Perl script. To get it up and running on your system there a few simple steps:

* Download the code!
* Make the main script accessible by adding it to your `PATH` with environment modules or other methods _(optional)_
* Set up the configuration file to work with your cluster environment
* Add some genome paths

## Downloading Cluster Flow
To get started with Cluster Flow, first download the source files using the links on the left. Stable releases are tagged as downloads. The main github branch will typically contain the latest stable release, with development work going on within development branches.

## Configuring your `PATH`
Linux systems use an environment variable called `PATH` to keep the locations of executable binaries. If you add the Cluster Flow directory to this, you'll be able to call the Cluster Flow `cf` file from anywhere.

### Environment modules
Many compute clusters use [environment modules](http://modules.sourceforge.net/) to manage the user’s `PATH`.

If your system uses environment modules, you can create a new module file to add Cluster Flow to the list of available modules. A typical module file example is given below - note that you'll need to change the file paths to where the Cluster Flow files are kept. The code below takes the version of Cluster Flow from the module script filename.

	#%Module1.0#####################################################################
	##
	## clusterflow modulefile
	##

	source $env(MODULE_INCLUDE)/functions.tcl
	getCluster

	set components [ file split [ module-info name ] ]
	set version [ lindex $components 1 ]

	set     modroot          /path/to/clusterflow/$version

	proc ModulesHelp { } {
	        global version modroot

	        puts stderr "\tclusterflow - use clusterflow $version"
	        puts stderr "\n\tVersion $version\n"
	}

	module-whatis   "Loads clusterflow environment."

	# Only one version at a time
	conflict clusterflow

	#Log loading to syslog
	logToSyslog

	if [module-info mode load] {
	    prepend-path        PATH            $modroot
	}

	if [module-info mode remove] {
	    remove-path         PATH            $modroot
	}

Once done, you may need to rebuild your module cache. You can load the new module by running

	module load clusterflow

You may also wish to create a symlinked directory called cf so that `module load cf` also works.

### Manually adding Cluster Flow to the `PATH`
You can manually add the directory containing cluster flow to your `PATH` with the following command :

	PATH=/path/to/clusterflow:$PATH

### Checking it's worked
For both of the above methods, you can check that the system can find Cluster Flow with the following command:

	which cf

You should see the directory containing the Cluster Flow files.

## Configuring Cluster Flow
### Cluster Flow config files
Cluster flow will search three locations for a config file every time it is run. Variables found in each file can override those read from a previous config file. They are, in order of priority:

* `<working directory>/clusterflow.config`
	* A config file found in the current working directory when a pipeline is executed has top priority, trumped only by command line parameters
* ``~/clusterflow.config``
	* A config file in your home directory can be used to set parameters such as notification level and e-mail address
* `<installation directory>/clusterflow.config`
	* A config file in the Cluster Flow installation directory is ideal for common settings specific to the environment

Config files contain key: value pairs. Syntax is as follows: `@key value` (tab delimited, one per line). Cluster Flow ships with an example config file called [clusterflow.config.example](https://github.com/ewels/clusterflow/blob/master/clusterflow.config.example "Browse the example config file on GitHub")

Typically, there will be a config file in the installation directory which contains the settings that make Cluster Flow work, then each user will have a personal configuration file in their home directory containing settings such as a notification e-mail address.

### Making Cluster Flow work with your environment
The key things to set up when installing Cluster Flow are the variables that dictate how CF should interact with your cluster - what commands it should use to submit jobs.

Cluster Flow currently supports GRIDEngine (SGE), SLURM and LSF. You can choose which one you're running with the following:

	/* Options: GRIDEngine, SLURM or LSF */
	@cluster_environment    SLURM

In most cases, that should be enough to get Cluster Flow to work! However, some people have some specific variables that need to be submitted with batch jobs (eg. project identifiers, time limits, other custom flags). If this is the case, the job submission command can be customised with the `@custom_job_submit_command` config variable.

To use this, enter your typical submission command with the following placeholders which will be replaced at run time:

* `{{ "{{command" }}}}`
	* The actual command which will be run to execute the module file
* `{{ "{{job_id" }}}}`
	* The unique identifier which will be assigned to use job dependencies
* `{{ "{{outfn" }}}}`
	* The filename of the log file to capture `STDOUT`
* `{{ "{{cores" }}}}`
	* How many cores to assign
* `{{ "{{mem" }}}}`
	* How much memory to assign
* `{{ "{{priority" }}}}`
	* A priority to set, defined in the config file
* `{{ "{{email" }}}}`
	* The user's e-mail address for cluster job notifications (if set)
* `{{ "{{notifications" }}}}`
	* A string describing which notifications to be sent (syntax depends on environment set above)

For example:
	
	@custom_job_submit_command      sbatch  -A MY_PROJECT_ID -t 2-00:00:00 -p core -n {{ "{{cores" }}}} --open-mode=append -o {{ "{{outfn" }}}} -J {{ "{{job_id" }}}} {{ "{{notifications" }}}} --wrap="{{ "{{command" }}}}"

Cluster Flow will generate it's own sensible default if this isn't set, so it's worth trying it without first.

**Note** - Cluster Flow will append the job dependency strings to the end of your custom command which are system specific, so it's important that `@cluster_environment` is still correct.

### System specific Gotchas
There are a few non-obvious things that have to be set correctly for some environments to work..

* GridEngine / SGE
  * GridEngine defaults to using `csh` to execute commands instead of the more standard `bash`. This will break some things. To override this behaviour, include the parameter `-S /bin/bash` for the commands


## Adding genome paths
Many modules within Cluster Flow require reference genomes for alignment. You can do this manually (probably quicker if you have lots to add) or you can use the interactive wizard that comes with Cluster Flow - just run:

	 cf --add_genome

Genome paths are stored within files called `genomes.config` which are stored in the same directories as `clusterflow.config`. Within this file, each path is described with `@genome_path` followed by a unique IDy used when submitting the Cluster Flow run (eg. `--genome GRCh37`). These are then followed by an absolute path. Optionally, species and assembly can be added after this (see [example](https://github.com/ewels/clusterflow/blob/master/genomes.config.example "Browse the example genomes file on GitHub")).

There are four types of paths that can be specified:

* `@genome_path`
	* The directory containing a reference genome
* `@bowtie_path`
	* The file name base of a set of bowtie indices
* `@bowtie2_path`
	* The file name base of a set of bowtie 2 indices
* `@gtf_path`
	* The file name of a GTF file for a given genome.

All four types of path should share genome keys if applicable. The fields should be separated by a tab character. Cluster Flow ships with an example genomes file called [genomes.config.example](https://github.com/ewels/clusterflow/blob/master/genomes.config.example "Browse the example genomes file on GitHub")

## Environment Module Aliases
If using environment modules, you may get some errors claiming that certain
tools are not installed. If you think that you do have that tool installed,
it could be because of a minor difference in the module name (eg. `fastqc`
versus `FastQC`). To avoid having to change the name of your modules, you
can configure aliases in your configuration file. You can also use these
aliases to specify specific software versions for Cluster Flow.

Aliases are added with the `@environment_module_alias` tag. For example:

	@environment_module_alias	fastqc	FastQC/0.11.2
	@environment_module_alias	trim_galore	TrimGalore

## Running Cluster Flow
Ok, you're done! Cluster Flow should now work (see [General Usage]({{site.baseurl}}/usage/) for instructions on how to run Cluster Flow). If you get any weird errors, have a look through the  [Troubleshooting]({{site.baseurl}}/troubleshooting/) section and if in doubt, drop me an e-mail. Good luck!

-------

# Config File reference
The following section describes the available variables that can be set in the config file. For an example, see the [clusterflow.config.example](https://github.com/ewels/clusterflow/blob/master/clusterflow.config.example "Browse the example config file on GitHub") file that comes bundled with Cluster Flow.

### @email
Sets your e-mail address, used for e-mail notifications.

### @check_updates
Cluster Flow can automatically check for new versions. If an update is available, it will print a notification each time you run a job. You can specify how often Cluster Flow should check for updates with this parameter. The syntax is a number followed by `d`, `w`, `m` or `y` for days, weeks, months or years. Cluster Flow will check for an update at runtime if this period or more has elapsed since you last ran it. You can disable update checks and alerts by setting `@check_updates 0` in your `~/clusterflow.config` file.

You can get Cluster Flow to manually check for updates by running cf `--check_updates`

### @split_files
The default number of input files to send to each run. Typically set to 1.

### @priority
The priority to give to cluster jobs.

### @max_runs
The maximum number of parallel runs that cluster flow will set off in one go. Default is 12 to avoid swamping the cluster for all other users.

### @notification
Multiple `@notification` key pairs can be set with the following values:

* complete
	* An e-mail notification is sent when all processing for all files has finished
* run
	* An e-mail is sent when each run finishes (each set of input files)
* end
	* A qsub notification e-mail is sent when each cluster job ends. Likely to result in a full inbox!
* suspend
	* A qsub notification e-mail is sent if a job is suspended
* abort
	* A qsub notification e-mail is sent if a job is aborted

Cluster Flow sends the `run` and `complete` notifications using the `cf_run_finished` and `cf_runs_all_finished` modules. These modules handle several tasks, such as cleaning useless warning messages from log files. E-mails contain the contents of all log files, plus a section at the top of highlighted messages, specified within log messages by being prefixed with `###CF`.

### @total_cores
The total number of cores available to a Cluster Flow pipeline. Modules are given a recommended number of cores so that resources can be allocated without swamping the cluster.

### @total_mem
The total amount of memory available to a Cluster Flow pipeline. Modules are given a recommended quota so that resources can be allocated without swamping the cluster.

### @cluster_environment
Cluster Flow can submit jobs to both GRIDEngine and LSF cluster management systems. This configuration variable sets which environment to use. The possible options are:

* GRIDEngine
* SLURM
* LSF

It is worth noting that Cluster Flow has been developed on GRIDEngine and SLURM Clusters so is likely to be most robust on similar setups. The `--qstat` and `--qstatall` parameters do not currently work with LSR. If you'd like to have a go at implementing it, that would be great!

### @ignore_modules
If you do not use environment modules on your system, you can prevent Cluster Flow from trying to use them (and giving a warning) by adding this line to your config file.

### @environment_module_alias
See above for instructions on how to use this tag.

### @log_highlight_string and @log_warning_string
To pull out specific highlights or warnings from log files, you can specify
search strings with these tags. If found, the e-mail will be highlighted
accordingly and the lines from the log file will be displayed at the top
of the report e-mail.

For example:

	@log_highlight_string at least one reported alignment
	@log_warning_string job failed
