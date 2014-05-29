---
title: Writing Pipelines & Modules
layout: toc
---

# Pipeline syntax
All pipelines conform to a standard syntax. The name of the pipeline is given by the filename, which should end in `.config`. The top of the file should contain a title and description surrounded by `/*` and `*/`

Variables are set using the same `@key value` syntax as in `clusterflow.config` files. Typical variables for pipelines are `@require_genome`, `@require_bowtie` or `@require_gtf`

Modules are described using `#` prefixes. Tab indentation denotes dependencies between modules. Syntax is `#module_name parameters`, where there can be any number of space separated parameters which will be passed on to the module at run time.

## Example pipeline
Here is an example pipeline, which requires a genome path and uses three modules:

	/*
	Example Pipeline
	================
	This pipeline is an example of running three modules which depend on each other. Module 2 is run with a parameter that modifies its behaviour. This block of text is used when cf --help example_pipeline is run
	*/
	#module1
	       #module2
	       #module2 parameter
	             #module3


Remember to run `dos2unix` on your pipeline before you run it, if you're working on a windows machine.

# Module syntax
Modules do not conform to any syntax _per se_ – they can be written in any language and in any style, though they must conform to a common API so that cluster flow can interact with them.

## Wrappers
If you have an existing script or tool, it's almost certain that it'll be easier to write a new Cluster Flow wrapper module than it will be to customise your script to be compatible with CF.

It's very easy to execute other code from within a wrapper perl script and will save you a lot of pain!

## Example module
An example module comes bunded with Cluster Flow, containing some typical pseudocode which you can modify for your own uses. You can view it on GitHub [here](https://github.com/ewels/clusterflow/blob/master/modules/example_module "Browse the Cluster Flow example module code on GitHub").

## Required command line flags
All modules will be called directly by Cluster Flow before any cluster jobs are set up, to determine their required job parameters. Each module must return a value for the following command line flags:

Flag | Description
-----|-------------
`--cores <num>` | Print required number of cores
`--mem <num>` | Print required amount of memory Print Help
`--modules` | Print names of required environment modules
`--help` | Print Help
{: .table}

### `--cores`
Cluster Flow will calculate how many cores are available for the module before launching it as a job on the cluster. It will call the module directly with the `--cores` flag and the suggested number of cores that the module can use. The module should print a single integer representing the number of cores that it requires and exit.

The maximum number of cores is a recommendation only. Cluster Flow will assign whatever the module returns.

### `--mem`
Much like cores, Cluster Flow will calculate the maximum memory available for the module and call the module with the `--mem` flag followed by an amount of memory in bytes. The module should print the required memory and exit. Required memory can be returned in bytes or formatted with a number followed by a single letter suffix; eg. `4G`, `4096M` or `4294967296`.

You can use the core Cluster Flow helper functions `allocate_cores()` and `allocate_memory()` to convert simple human readable strings to bytes. For example:

	return CF::Helpers:: allocate_memory($allocated, '2G', '10G');

The maximum available memory is a recommendation only. Cluster Flow will assign whatever the module returns.


### `--modules`
Most Cluster Flow modules will require a system program to exist in the `PATH`. It’s common practice to use environment modules to manage this. Environment modules need to be loaded from the head node before the qsub jobs are set off. Cluster Flow calls each CF module with the `--modules` flag before running it. A comma separated list of module names should be printed. Each of these will be loaded with `module load <name>`.
	
### `--help`
Cluster Flow can be called with the `--help` flag followed by any pipeline or module name. If a module name is requested, Cluster Flow will call that module with the `--help` parameter and print the `STDOUT`.

## Command line parameters
When a module is called by its cluster job, it will have the following parameters passed in this order:

* Run file filename
* Cluster job ID
* Previous job ID
* Number of cores assigned
* Amount of memory assigned
* Extra parameters

A run file is created by Cluster Flow for each batch of files. It describes variables to be used, the pipeline specified and the filenames used by each module. The syntax of variables and pipeline is described in Pipeline syntax.

File names are described by a job identifier followed by a tab then a filename. Each module is provided with its own job ID and the ID of the job that was run previously. By using these identifiers, the module can read which input files to use and write out the resulting filenames to the run file when complete. Example run file syntax:

	first_job_938 filename_1.txt
	first_job_938 filename_2.txt
	second_job_375 filename_1_processed.txt
	second_job_375 filename_2_processed.txt

There can be any number of extra parameters, these are specific the module and are specified in the pipeline configuration.

## E-mail report highlights
Any `STDOUT` or `STDERR` that your module produces will be written to the Cluster Flow log file. At the end of each run and pipeline, an e-mail will be sent to the submitter with details of the run results (if specified by the config settings). Because the log file can be very long Cluster Flow pulls out any lines starting with `###CF`. Typically, such a line should be printed when a module finishes, with a concise summary of whether it worked or not.

Modules should print the command that they are going to run to STDERR so that this is recorded in the log file. These are also sent in the e-mail notification and should start with `###CFCMD`

## Exit codes
It's likely that your cluster will continue to fire off the dependent jobs as soon as the parent jobs finish, irrespective of their output. If a module fails, the cleanest way to exit is with a success code, but without printing any resulting output filename. The following modules will not find their input filenames and so should immediately exit.

## Example module
An example module is distributed with Cluster Flow in `modules/example_module`

The module doesn’t do anything but shows the typical workflow for a Perl module, replete with lots of comments saying what everything is doing. Typically, the easiest way to write a new module is to copy and old one, and modify it to suit your needs.

If you write a new module, please let us know! We’d love to package it with Cluster Flow so that other people can benefit.

# Helper Perl Packages
If your module is written in Perl, there are some common Cluster Flow packages that you can use to provide some pre-written functions.

Ok, so the heading should be "Perl Modules", but I thought that would be too confusing. I'll stick with packages from now on too.

## Using Cluster Flow packages
There are currently three packages available to Cluster Flow modules. `Helpers` contains subroutines of general use for most modules. `Constants` and `Headnodehelpers` contain subroutines primarily for use in the main `cf` script. You can include the Helpers package by adding the following to the top of your module file:

	use FindBin qw($Bin);
	use lib "$FindBin::Bin/../source";
	use CF::Helpers;

We use the package `FindBin` to add the binary directory to the path (where `cf` is executing from).

## Cluster Flow helper functions
### `load_runfile_params(@ARGV);`
This function reads the run file and returns a list of input filenames. Note that `$files` is an array reference and `$config` is a hash reference.

	my ($files, $runfile, $job_id, $prev_job_id, $cores, $mem, $parameters, $config) = CF::Helpers::load_runfile_params(@ARGV);

### `is_paired_end(@files);`
This function takes an array of file names and returns an array of single end files and an array of arrays of paired end files. It sorts the files alphabetically, removes any occurance of `_[1-4]` from the filename and compares. Identical pairs are returned as paired end. This function will return the appropriate file arrays if `--paired` or `--single` was used at run time.

	my ($se_files, $pe_files) = CF::Helpers::is_paired_end(@$files); foreach my $file (@$se_files){
	       print “$file is single end\n”;
	}
	foreach my $files_ref (@$pe_files){
		my @files = @$files_ref;
		print $files[0].” and “.$files[1].” are paired end.\n”;
	}

### `is_bam_paired_end($files);`
Looks at BAM/SAM file headers and tries to determine whether it has been generated using paired end input files or single end. The subroutine simply searches the `@PG` header for the occurrence of `-1` and `-2`.

	if(CF::Helpers::is_bam_paired_end($file)){
		# do something with paired end BAM
	} else {
		# do something with single end BAM
	}

### `fastq_encoding_type($file);`
Scans a FastQ file and tries to determine the encoding type. Returns strings `integer`, `solexa`, `phred33`, `phred64` or `0` if too few reads to safely determine. This is done by observing the minimum and maximum quality scores.

For more details, see the [Wikipedia page on FastQ encoding](http://en.wikipedia.org/wiki/FASTQ_format#Encoding)

	($encoding) = CF::Helpers::fastq_encoding_type($file);

### `fastq_min_length($file);`
Scans the first 100000 reads of a FastQ file and returns the longest read length that it finds.

### `parse_seconds($raw, [$long]);`
Simple function that takes time in seconds as an input and returns a human readable string. The optional second `$long` variable determines whether to use h/m/s (`0`, false) or hours/minutes/seconds (`1`, true) and defaults to true.

### `human_readable_to_bytes($raw);`
### `bytes_to_human_readable($raw);`
Two functions which convert between human readable memory strings (eg. `4G` or `100M`) and bytes.

### `allocate_cores($recommended, $min, $max);`
Takes the suggested number of cores to use, a minimum and maximum number and returns a sensible result.

### `allocate_memory($recommended, $min, $max);`
Takes the suggested number of memory to use, a minimum and maximum amount and returns a sensible result. Input can be human readable strings or bytes. Returns a value in bytes.

