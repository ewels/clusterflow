## Cluster Flow config files
Cluster flow will search three locations for a config file every time it is run. Variables found in each file can override those read from a previous config file. They are, in order of priority:

* `<working directory>/clusterflow.config`
	* A config file found in the current working directory when a pipeline is executed has top priority, trumped only by command line parameters
* ``~/clusterflow.config``
	* A config file in your home directory can be used to set parameters such as notification level and e-mail address
* `<installation directory>/clusterflow.config`
	* A config file in the Cluster Flow installation directory is ideal for common settings specific to the environment

Config files contain key: value pairs. Syntax is as follows: `@key value` (tab delimited, one per line). Cluster Flow ships with an example config file called `clusterflow.config.example`

Typically, there will be a config file in the installation directory which contains the settings that make Cluster Flow work, then each user will have a personal configuration file in their home directory containing settings such as a notification e-mail address.

### Making Cluster Flow work with your environment
The key things to set up when installing Cluster Flow are the variables that dictate how CF should interact with your cluster - what commands it should use to submit jobs.

Cluster Flow currently supports GRIDEngine (SGE), SLURM and LSF. You can choose which one you're running with the following:

	/* Options: GRIDEngine, SLURM or LSF */
	@cluster_environment    SLURM

In most cases, that should be enough to get Cluster Flow to work! However, some people have some specific variables that need to be submitted with batch jobs (eg. project identifiers, time limits, other custom flags). If this is the case, the job submission command can be customised with the `@custom_job_submit_command` config variable.

To use this, enter your typical submission command with the following placeholders which will be replaced at run time:

* `{{command}}`
	* The actual command which will be run to execute the module file
* `{{job_id}}`
	* The unique identifier which will be assigned to use job dependencies
* `{{outfn}}`
	* The filename of the log file to capture `STDOUT`
* `{{cores}}`
	* How many cores to assign
* `{{mem}}`
	* How much memory to assign
* `{{priority}}`
	* A priority to set, defined in the config file
* `{{email}}`
	* The user's e-mail address for cluster job notifications (if set)
* `{{notifications}}`
	* A string describing which notifications to be sent (syntax depends on environment set above)

For example:

	@custom_job_submit_command      sbatch  -A MY_PROJECT_ID -t 2-00:00:00 -p core -n {{cores}} --open-mode=append -o {{outfn}} -J {{job_id}} {{notifications}} --wrap="{{command}}"

Cluster Flow will generate it's own sensible default if this isn't set, so it's worth trying it without first.

**Note** - Cluster Flow will append the job dependency strings to the end of your custom command which are system specific, so it's important that `@cluster_environment` is still correct.

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

### Running Cluster Flow
Ok, you're done! Cluster Flow should now work (see [General Usage](usage/) for instructions on how to run Cluster Flow). If you get any weird errors, have a look through the  [Troubleshooting](troubleshooting/) section and if in doubt, drop me an e-mail. Good luck!

-------

## Config File reference
The following section describes the available variables that can be set in the config file. For an example, see the `clusterflow.config.example` file that comes bundled with Cluster Flow.

#### @email
Sets your e-mail address, used for e-mail notifications.

#### @check_updates
Cluster Flow can automatically check for new versions. If an update is available, it will print a notification each time you run a job. You can specify how often Cluster Flow should check for updates with this parameter. The syntax is a number followed by `d`, `w`, `m` or `y` for days, weeks, months or years. Cluster Flow will check for an update at runtime if this period or more has elapsed since you last ran it. You can disable update checks and alerts by setting `@check_updates 0` in your `~/clusterflow.config` file.

You can get Cluster Flow to manually check for updates by running cf `--check_updates`

#### @split_files
The default number of input files to send to each run. Typically set to 1.

#### @priority
The priority to give to cluster jobs.

#### @max_runs
The maximum number of parallel runs that cluster flow will set off in one go. Default is 12 to avoid swamping the cluster for all other users.

#### @notification
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

#### @total_cores
The total number of cores available to a Cluster Flow pipeline. Modules are given a recommended number of cores so that resources can be allocated without swamping the cluster.

#### @total_mem
The total amount of memory available to a Cluster Flow pipeline. Modules are given a recommended quota so that resources can be allocated without swamping the cluster.

#### @cluster_environment
Cluster Flow can submit jobs to both GRIDEngine and LSF cluster management systems. This configuration variable sets which environment to use. The possible options are:

* GRIDEngine
* SLURM
* LSF

It is worth noting that Cluster Flow has been developed on GRIDEngine and SLURM Clusters so is likely to be most robust on similar setups. The `--qstat` and `--qstatall` parameters do not currently work with LSR. If you'd like to have a go at implementing it, that would be great!

#### @ignore_modules
If you do not use environment modules on your system, you can prevent Cluster Flow from trying to use them (and giving a warning) by adding this line to your config file.

#### @environment_module_alias
See above for instructions on how to use this tag.

#### @log_highlight_string and @log_warning_string
To pull out specific highlights or warnings from log files, you can specify
search strings with these tags. If found, the e-mail will be highlighted
accordingly and the lines from the log file will be displayed at the top
of the report e-mail.

For example:

	@log_highlight_string at least one reported alignment
	@log_warning_string job failed


-----



## Command Line Reference


Flag                      | Description
------------------------- | ---------------------------------------------------
`--genome <ID>`           | ID of a genome referred to in `genomes.config `
`--genome_path <path>`    | Path to a genome to be used for alignment
`--bowtie_path <path>`    | Path to a bowtie index base to be used for alignment
`--gtf_path <path>`       | Path to a GTF file to be used for alignment (eg. for Tophat)
`--paired`                | Force paired-end mode
`--single`                | Force single-end mode
`--no_fn_check`           | Disable input file type checking
`--file_list`             | Text file containing input files or download URLs
`--params`                | Specify extra module parameters for this run
`--split_files <num>`     | Create one run per `<num>` files
`--max_runs <num>`        | Divide input files into `<num>` runs. Set as 0 to disable.
`--email <email>`         | Set the e-mail address for notifications
`--priority <num>`        | Set the queue priority for cluster jobs
`--cores <num>`           | Set the maximum number of cores to use for all runs
`--mem <string>`          | Set the maximum memory to use for all runs
`--notifications [cresa]` | Specify desired notifications
`--list_pipelines`        | Print available pipelines
`--list_modules`          | Print available modules
`--list_genomes`          | Print available genomes
`--dry_run`               | Prints jobs to terminal instead of submitting them to the cluster
`--qstat`                 | Displays formatted qstat output of your jobs
`--qstatall`              | Displays formatted qstat output of all jobs  
`--qstatcols`             | Colours output from --qstat or --qstatall
`--qdel <id>`             | Delete all jobs from a running pipeline. `<id>` is printed with --qstat
`--make_config`           | Interactive prompt to generate a personalised CF config file
`--add_genome`            | Interactive wizard to add new genomes to your genomes.config files
`--version`               | Print version of Cluster Flow installed
`--check_updates`         | Look for available Cluster Flow updates
`--help`                  | Print help message

## Parameter Details
### `--genome`
**Default: none**

Some pipelines which carry out a reference genome alignment require a genome directory path to be set. Requirements for format may vary between modules.

### `--paired`
**Default: Auto-detect**

If specified, Cluster Flow will send two files to each run, assuming that the order that the file list is supplied in corresponds to two read files. If an odd number of files is supplied, the final file is submitted as single end.

### `--single`
**Default: Auto-detect**

If specified, Cluster Flow will ignore its auto-detection of paired end input files and force the single end processing of each input file.

### `--no_fn_check`
**Default: none**

Cluster Flow will make sure that all of the input files have the same file extension to avoid accidentally submitting files that aren’t part of the run. Specifying this parameter disables this check.

### `--file_list`
**Default: none**

If specified, you can define a file containing a list of filenames to pass to the pipeline (one per line). This is particularly useful when supplying a list of download URLs.

### `--params`
**Default: none**

Pipelines and their modules are configured to run with sensible defaults. Some modules accept parameters which change their behaviour. Typically, these are set within a pipeline config file. By using `--params`, you can add extra parameters at run time. These will be set for every module in the pipeline (though they probably won’t all recognise them).

### `--split_files`
**Default: (config file - typically 1)**

Cluster Flow generates multiple parallel runs for the supplied input files when run. This is typically a good thing, the cluster is designed to run jobs in parallel. Some jobs may involve many small tasks with a large number of input files however, and 1:1 parallelisation may not be practical. In such cases, the number of input files to assign to each run can be set this flag.

### `--max_runs`
**Default: none**

It can sometimes be a pain to count the number of input files and work out a sensible number to use with `--split_files`. Cluster Flow can take the `--max_runs` value and divide the input files into this number of runs, setting `--split_files` automatically.

A default can be set for `--max_runs` in the `clusterflow.config` file, and this value is set to 12 if no value is found in the config files. Set to 0 to disable.

This parameter will override anything set using `--split_files`.

### `--email`
**Default: (config file)**

Cluster Flow can send notification e-mails regarding the status of runs. Typically, e-mail
address should be set using a personalised ~/clusterflow.config value (see below). This parameter allows you to override that setting on a one-off basis.

### `--priority`
**Default: (config file - typically -500)**

Many cluster managers can use a priority system to manage jobs in the queue. Typically, GRIDEngine priorities can be set ranging from -1000 to 0.

### `--cores`
**Default: (config file - typically 64)**

Override the maximum number of cores allowed for each Cluster Flow pipeline, typically set in the Cluster Flow config file. For more information see Avoiding cluster overload.

### `--mem`
**Default: (config file - typically 128G)**

Setting `--mem` allows you to override the maximum amount of simultaneously assigned memory. For more information see Avoiding cluster overload.

### `--notifications`
**Default: (config file - typically cea)**

Cluster Flow can e-mail you notifications about the progress of your runs. There are several levels of notification that you can choose using this flag. They are:

* `c` - Send notification when all runs in a pipeline are completed
* `r` - Send a notification when each run is completed
* `e` - Send a notification when a cluster job ends
* `s` - Send a notification if a cluster job is suspended
* `a` - Send a notification if a cluster job is aborted

Setting these options at run time with the `--notifications` flag will override the settings present in your clusterflow.config configuration files.
Note: setting the `s` flag when using many input files with a long pipeline may cause your inbox to be flooded.

### `--qstat` and `--qstatall`
When you have a lot of jobs running and queued, the qstat summary can get a little overwhelming. To combat this and show job hierarchy in an intuitive manner, you can enter into the console `cf --qstat`. This parses qstat output and displays it nicely. `cf --qstatall` does the same but for all jobs by all users.

You'll probably find that you want to run this command quite a lot. To make it a little less clumsy, you can create aliases in your `.bashrc` script. Typically, these two lines:

	alias qs='cf --qstat'
	alias qsa='cf --qstatall'

If you’re feeling lazy, you can append these lines to your .bashrc script through the command line by copying and pasting the following commands:

```
sed -i "$ a\alias qs='cf --qstat'" ~/.bashrc
sed -i "$ a\alias qsa='cf --qstatall'" ~/.bashrc
```

Note - these guys don't work with LSF yet, writing the code to parse the status command outputs is a bit of a pain when I don't have a LSF testing server to work on. If anyone fancies contributing some code to do so, that would be great!

### `--qstatcols`
How this will look depends entirely on your terminal colour setup. It works nicely in terminal windows with light backgrounds and can look really horrible on terminal windows with dark backgrounds. For this reason it is only enabled when specified. The author uses the lovely  [Solarized](http://ethanschoonover.com/solarized) theme by Ethan Schoonover and it looks nice when he runs it. Maybe you should use Solarized too ;)

If you like how it looks with colour, you can add the `--qstatcols` flag to the `.bashrc` aliases above to use them every time.

### `--qdel`
Sometimes you may be running multiple pipelines and want to stop just one. It can be a pain to find the job numbers to do this manually, so instead you can use Cluster Flow to kill these jobs. When running `cf --qstat`, ID values are printed for each pipeline. Use this with `--qdel`. eg:

```
cf --qdel sra_bowtie_1391074179
```
