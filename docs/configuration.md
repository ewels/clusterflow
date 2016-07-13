## Config file locations
Cluster flow will search three locations for a config file every time it is
run. Variables found in each file can override those read from a previous
config file. They are, in order of priority:

* `<working directory>/clusterflow.config`
	* A config file found in the current working directory when a pipeline is
      executed has top priority, trumped only by command line parameters.
* `~/clusterflow.config`
	* A config file in your home directory can be used to set parameters such
      as notification level and e-mail address.
* `<installation directory>/clusterflow.config`
	* A config file in the Cluster Flow installation directory is ideal for
      common settings specific to the environment.

Config files contain key: value pairs. Syntax is as follows: `@key value`
(tab delimited, one per line). The Cluster Flow source code comes with an
example config file called
[`clusterflow.config.example`](https://github.com/ewels/clusterflow/blob/master/clusterflow.config.example)

Typically, there will be a config file in the installation directory which
contains the settings that make Cluster Flow work. Each user will then have a
personal configuration file in their home directory containing settings such
as a notification e-mail address.

## Environment Setup
The key things to set up when installing Cluster Flow are the variables that
dictate how CF should interact with your cluster - what commands it should use
to submit jobs.

Cluster Flow currently supports GRIDEngine (SGE), SLURM and LSF, as well as 
running locally using background bash jobs. You can specify which environment
to use with `@cluster_environment`:
```
/* Options: local, GRIDEngine, SLURM or LSF */
@cluster_environment    SLURM
```

In most cases, that should be enough to get Cluster Flow to work! However,
some people have some specific variables that need to be submitted with
batch jobs (eg. project identifiers, time limits, other custom flags). If this
is the case, the job submission command can be customised with the
`@custom_job_submit_command` config variable.

To use this, enter your typical submission command with the following
placeholders which will be replaced at run time:

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
* `{{time}}`
	* How much time to assign
* `{{priority}}`
	* A priority to set, defined in the config file
* `{{project}}`
	* The cluster project to use
* `{{qname}}`
    * The cluster queue name
* `{{email}}`
	* The user's e-mail address for cluster job notifications (if set)
* `{{notifications}}`
	* A string describing which notifications to be sent (syntax depends on
      environment set above)

Simply omit any variables which are not needed on your cluster. For example:
```
@custom_job_submit_command      sbatch  -A MY_PROJECT_ID -t 2-00:00:00 -p core -n {{cores}} --open-mode=append -o {{outfn}} -J {{job_id}} {{notifications}} --wrap="{{command}}"
```

Cluster Flow will generate it's own sensible default if this isn't set, so
it's worth trying it without first.

> _Note:_ Cluster Flow will append the job dependency strings to the end of
> your custom command which are system specific, so it's important that
> `@cluster_environment` is correct.

## Config File reference
The following section describes the available variables that can be set in the
config file. For an example, see the `clusterflow.config.example` file that
comes bundled with Cluster Flow.

### @email
Sets your e-mail address, used for e-mail notifications.

### @colourful / @colorful
Set to true to make the output from `cf --qstat` and `cf --qstatall`
colourful (and hopefully easier to read).
```
@colourful	1
```

### @merge_regex
A regex used to automatically merge files before pipeline processing starts.
This works by matching a single regex group within a filename. If multiple
input files have the same matching group, they will be merged. The regex
group is then used to give the output filename.

For example, given the following config regex:
```
@merge_regex    [1-8]_[0-9]{6}_[a-zA-Z0-9]+_(P\d+_\d+_[12]).fastq.gz
```

These input files:
```
1_160312_CDSH32SDB3889_P1234_001_1.fastq.gz
1_160312_CDSH32SDB3889_P1234_001_2.fastq.gz
2_160312_CDSH32SDB3889_P1234_001_1.fastq.gz
2_160312_CDSH32SDB3889_P1234_001_2.fastq.gz
```
Would give the resulting merged files:
```
P1234_001_1.fastq.gz
P1234_001_2.fastq.gz
```

### @split_files
The default number of input files to send to each run. Typically set to 1.

### @max_runs
The maximum number of parallel runs that cluster flow will set off in one go. Default is 12 to avoid swamping the cluster for all other users.

### @total_cores
The total number of cores available to a Cluster Flow pipeline. Modules are
given a recommended number of cores so that resources can be allocated without
swamping the cluster.

### @total_mem
The total amount of memory available to a Cluster Flow pipeline. Modules are
given a recommended quota so that resources can be allocated without swamping
the cluster.

### @max_time
The maximum time that a job should request in a Cluster Flow pipeline. For
example, to prevent jobs from requesting more than 10 days:
```
@max_time	10-00
```

### @time_multiplier
If your cluster is running slowly and the default time limits specified
in Cluster Flow modules are not enough, jobs will fail due to timing out.
`@time_multiplier` is a quick and dirty way to avoid this. Setting 
`@time_multiplier` to `2` will double the requested time for every job.
Note that these times will still be capped by `@max_time`.

### @priority
The priority to give to cluster jobs.

### @cluster_environment / @custom_job_submit_command
See above docs: _[Environment setup](#environment-setup)_.

### @ignore_modules
If you do not use environment modules on your system, you can prevent Cluster
Flow from trying to use them (and giving a warning) by adding this line to
your config file.

### @environment_module_always
Specify an environment module to always load for every Cluster Flow pipeline.
Can be used multiple times.

### @environment_module_alias
If using environment modules, you may get some errors claiming that certain
tools are not installed. If you think that you do have that tool installed,
it could be because of a minor difference in the module name (eg. `fastqc`
versus `FastQC`). You can configure aliases in your configuration file.
You can also use these aliases to specify specific software versions for
Cluster Flow.

Aliases are added with the `@environment_module_alias` tag. For example:
```
@environment_module_alias	fastqc	FastQC/0.11.2
@environment_module_alias	trim_galore	TrimGalore
```

### @log_highlight_string / @log_warning_string
To pull out specific highlights or warnings from log files, you can specify
search strings with these tags. If found, the e-mail will be highlighted
accordingly and the lines from the log file will be displayed at the top
of the report e-mail.

For example:
```
@log_highlight_string at least one reported alignment
@log_warning_string job failed
```

### @notification
Multiple `@notification` key pairs can be set with the following values:

* `complete`
	* A Cluster Flow e-mail notification is sent when all processing for all
      files has finished
* `run`
	* A Cluster Flow e-mail is sent when each run finishes (each set of input files)
* `end`
	* A cluster notification e-mail is sent when each cluster job ends. Likely
      to result in a full inbox!
* `suspend`
	* A cluster notification e-mail is sent if a job is suspended
* `abort`
	* A cluster notification e-mail is sent if a job is aborted

Cluster Flow sends the `run` and `complete` notifications using the
`cf_run_finished` and `cf_runs_all_finished` modules. These modules handle
several tasks, such as cleaning useless warning messages from log files.
E-mails contain the contents of all log files, plus a section at the top
of highlighted messages, specified within log messages by being prefixed
with `###CF`.

### @check_updates
Cluster Flow can automatically check for new versions. If an update is
available, it will print a notification each time you run a job. You can
specify how often Cluster Flow should check for updates with this parameter.
The syntax is a number followed by `d`, `w`, `m` or `y` for days, weeks,
months or years. Cluster Flow will check for an update at runtime if this
period or more has elapsed since you last ran it. You can disable update
checks and alerts by setting `@check_updates 0` in your
`~/clusterflow.config` file.

You can manually get Cluster Flow to check for updates by running
`cf --check_updates`
