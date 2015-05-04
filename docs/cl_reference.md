---
title: Command Line Reference
layout: toc
---

# Command Line Parameter Overview


Flag | Description
---- | -----------
`--genome <ID>` | ID of a genome referred to in `genomes.config `
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
`--qdel <id>` | Delete all jobs from a running pipeline. `<id>` is printed with --qstat
`--make_config` | Interactive prompt to generate a personalised CF config file 
`--add_genome` | Interactive wizard to add new genomes to your genomes.config files 
`--version` | Print version of Cluster Flow installed 
`--check_updates` | Look for available Cluster Flow updates 
`--help` | Print help message
{: .table}

# Parameter Details
## `--genome`
**Default: none**

Some pipelines which carry out a reference genome alignment require a genome directory path to be set. Requirements for format may vary between modules.

## `--paired`
**Default: Auto-detect**

If specified, Cluster Flow will send two files to each run, assuming that the order that the file list is supplied in corresponds to two read files. If an odd number of files is supplied, the final file is submitted as single end.

## `--single`
**Default: Auto-detect**

If specified, Cluster Flow will ignore its auto-detection of paired end input files and force the single end processing of each input file.

## `--no_fn_check`
**Default: none**

Cluster Flow will make sure that all of the input files have the same file extension to avoid accidentally submitting files that aren’t part of the run. Specifying this parameter disables this check.

## `--file_list`
**Default: none**

If specified, you can define a file containing a list of filenames to pass to the pipeline (one per line). This is particularly useful when supplying a list of download URLs.

## `--params`
**Default: none**

Pipelines and their modules are configured to run with sensible defaults. Some modules accept parameters which change their behaviour. Typically, these are set within a pipeline config file. By using `--params`, you can add extra parameters at run time. These will be set for every module in the pipeline (though they probably won’t all recognise them).

## `--split_files`
**Default: (config file - typically 1)**

Cluster Flow generates multiple parallel runs for the supplied input files when run. This is typically a good thing, the cluster is designed to run jobs in parallel. Some jobs may involve many small tasks with a large number of input files however, and 1:1 parallelisation may not be practical. In such cases, the number of input files to assign to each run can be set this flag.

## `--max_runs`
**Default: none**

It can sometimes be a pain to count the number of input files and work out a sensible number to use with `--split_files`. Cluster Flow can take the `--max_runs` value and divide the input files into this number of runs, setting `--split_files` automatically.

A default can be set for `--max_runs` in the `clusterflow.config` file, and this value is set to 12 if no value is found in the config files. Set to 0 to disable.

This parameter will override anything set using `--split_files`.

## `--email`
**Default: (config file)**

Cluster Flow can send notification e-mails regarding the status of runs. Typically, e-mail
address should be set using a personalised ~/clusterflow.config value (see below). This parameter allows you to override that setting on a one-off basis.

## `--priority`
**Default: (config file - typically -500)**

Many cluster managers can use a priority system to manage jobs in the queue. Typically, GRIDEngine priorities can be set ranging from -1000 to 0.

## `--cores`
**Default: (config file - typically 64)**

Override the maximum number of cores allowed for each Cluster Flow pipeline, typically set in the Cluster Flow config file. For more information see Avoiding cluster overload.

## `--mem`
**Default: (config file - typically 128G)**

Setting `--mem` allows you to override the maximum amount of simultaneously assigned memory. For more information see Avoiding cluster overload.

## `--notifications`
**Default: (config file - typically cea)**

Cluster Flow can e-mail you notifications about the progress of your runs. There are several levels of notification that you can choose using this flag. They are:

* `c` - Send notification when all runs in a pipeline are completed
* `r` - Send a notification when each run is completed
* `e` - Send a notification when a cluster job ends
* `s` - Send a notification if a cluster job is suspended
* `a` - Send a notification if a cluster job is aborted

Setting these options at run time with the `--notifications` flag will override the settings present in your clusterflow.config configuration files.
Note: setting the `s` flag when using many input files with a long pipeline may cause your inbox to be flooded.

## `--qstat` and `--qstatall` 
When you have a lot of jobs running and queued, the qstat summary can get a little overwhelming. To combat this and show job hierarchy in an intuitive manner, you can enter into the console `cf --qstat`. This parses qstat output and displays it nicely. `cf --qstatall` does the same but for all jobs by all users.

You'll probably find that you want to run this command quite a lot. To make it a little less clumsy, you can create aliases in your `.bashrc` script. Typically, these two lines:

	alias qs='cf --qstat'
	alias qsa='cf --qstatall'

If you’re feeling lazy, you can append these lines to your .bashrc script through the command line by copying and pasting the following commands:

	sed -i "$ a\alias qs='cf --qstat'" ~/.bashrc
	sed -i "$ a\alias qsa='cf --qstatall'" ~/.bashrc

Note - these guys don't work with LSF yet, writing the code to parse the status command outputs is a bit of a pain when I don't have a LSF testing server to work on. If anyone fancies contributing some code to do so, that would be great!

## `--qstatcols`
How this will look depends entirely on your terminal colour setup. It works nicely in terminal windows with light backgrounds and can look really horrible on terminal windows with dark backgrounds. For this reason it is only enabled when specified. The author uses the lovely  [Solarized](http://ethanschoonover.com/solarized) theme by Ethan Schoonover and it looks nice when he runs it. Maybe you should use Solarized too ;)

If you like how it looks with colour, you can add th `--qstatcols` flag to the `.bashrc` aliases above to use them every time.

## `--qdel`
Sometimes you may be running multiple pipelines and want to stop just one. It can be a pain to find the job numbers to do this manually, so instead you can use Cluster Flow to kill these jobs. When running `cf --qstat`, ID values are printed for each pipeline. Use this with `--qdel`. eg:

	cf --qdel sra_bowtie_1391074179