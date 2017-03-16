## Cluster Flow Command Line Reference
Cluster Flow pipelines are launched as follows:
```
cf [flags] <pipeline> <input-files>
```
These flags are used to customise run-time parameters for the pipeline that
Cluster Flow will launch.

### --genome
_Default: none_

Some pipelines which carry out a reference genome alignment require a genome
directory path to be set. Requirements for format may vary between modules.

### --paired
_Default: Auto-detect_

If specified, Cluster Flow will send two files to each run, assuming that the
order that the file list is supplied in corresponds to two read files. If an
odd number of files is supplied, the final file is submitted as single end.

### --single
_Default: Auto-detect_

If specified, Cluster Flow will ignore its auto-detection of paired end input
files and force the single end processing of each input file.

### --no_fn_check
_Default: none_

Cluster Flow will make sure that all of the input files have the same file
extension to avoid accidentally submitting files that aren’t part of the run.
Specifying this parameter disables this check.

### --file_list
_Default: none_

If specified, you can define a file containing a list of filenames to pass to
the pipeline (one per line). This is particularly useful when supplying a list
of download URLs.

### --params
_Default: none_

Pipelines and their modules are configured to run with sensible defaults. Some
modules accept parameters which change their behaviour. Typically, these are
set within a pipeline config file. By using `--params`, you can add extra
parameters at run time. These will be set for every module in the pipeline
(though they probably won’t all recognise them).

### --split_files
_Default: (config file - typically 1)_

Cluster Flow generates multiple parallel runs for the supplied input files when
run. This is typically a good thing, the cluster is designed to run jobs in
parallel. Some jobs may involve many small tasks with a large number of input
files however, and 1:1 parallelisation may not be practical. In such cases, the
number of input files to assign to each run can be set this flag.

### --max_runs
_Default: none_

It can sometimes be a pain to count the number of input files and work out a
sensible number to use with `--split_files`. Cluster Flow can take the
`--max_runs` value and divide the input files into this number of runs, setting
`--split_files` automatically.

A default can be set for `--max_runs` in the `clusterflow.config` file, and
this value is set to 12 if no value is found in the config files. Set to `0` to
disable.

This parameter will override anything set using `--split_files`.

### --runfile_prefix
_Default: none_

Optional custom prefix for run file filenames. This is useful if you are
running multiple instances of Cluster Flow with the same input file in
the same directory, as it avoids potential clashes / mixups. For example:
```
cf --runfile_prefix bt1 --genome GRCh37 fastq_bowtie1 my_sample.fq
cf --runfile_prefix bt2 --genome GRCh37 fastq_bowtie2 my_samplefq
```

### --ref
_Default: none_

Specify a reference genome without adding to the `genomes.config` file.
Should be in the format `<ref_type>=<path>`, eg:
```
cf --ref bowtie=/path/to/bowtie/index <pipeline> <files>
```

### --dry_run
_Default: false_

Do everything except for actually launching cluster jobs. Useful for
testing and checking that jobs will be created properly.

## Customising Behaviour
Typically, Cluster Flow settings are set in static configuration files.
However, sometimes it can be useful to specify parameters on a one-off
basis on the command line.

### --email
_Default: (config file)_

Cluster Flow can send notification e-mails regarding the status of runs.
Typically, e-mail address should be set using `@email` in `~/clusterflow.config`
(see above). This parameter allows you to override that setting on a one-off
basis.

### --priority
_Default: (config file - typically -500)_

Many cluster managers can use a priority system to manage jobs in the queue.
Typically, GRIDEngine priorities can be set ranging from -1000 to 0.

### --cores
_Default: (config file - typically 64)_

Override the maximum number of cores allowed for each Cluster Flow pipeline,
typically set in the Cluster Flow config file. For more information see
Avoiding cluster overload.

### --mem
_Default: (config file - typically 128G)_

Setting `--mem` allows you to override the maximum amount of simultaneously
assigned memory. For more information see Avoiding cluster overload.

### --time
_Default: (config file - typically none)_

Override the maximum requested time assigned to jobs. For more information,
see Avoiding cluster overload.

### --project
_Default: (config file - typically none)_

Specify the project to use on the cluster for this run.

### --qname
_Default: (config file - typically none)_

Specify a custom cluster queue to use for this run.

### --environment
_Default: (config file - custom)_

Override the default environment to use for this pipeline run. Useful for
testing or small jobs, can run using bash commands instead of submitting
cluster jobs. For example:
```
cf --environment local test_pipeline *.txt
```

### --notifications
_Default: (config file - typically cea)_

Cluster Flow can e-mail you notifications about the progress of your runs.
There are several levels of notification that you can choose using this flag.
They are:

* `c` - Send notification when all runs in a pipeline are completed
* `r` - Send a notification when each run is completed
* `e` - Send a notification when a cluster job ends
* `s` - Send a notification if a cluster job is suspended
* `a` - Send a notification if a cluster job is aborted

Setting these options at run time with the `--notifications` flag will override
the settings present in your clusterflow.config configuration files.
Note: setting the `s` flag when using many input files with a long pipeline may
cause your inbox to be flooded.

## Other Functions
These flags instruct Cluster Flow to do something other than submit a
pipeline.

### --qstat / --qstatall
When you have a lot of jobs running and queued, the qstat summary can get a
little overwhelming. To combat this and show job hierarchy in an intuitive
manner, you can enter into the console `cf --qstat`. This parses qstat output
and displays it nicely. `cf --qstatall` does the same but for all jobs by all
users.

You'll probably find that you want to run this command quite a lot. To make it
a little less clumsy, you can create aliases in your `.bashrc` or
`.bash_profile` scripts, which run every time you log in.
```
alias qs='cf --qstat'
alias qsa='cf --qstatall'
```

To append these lines to your `.bashrc` script you can use the following
command:
```
echo -e "alias qs='cf --qstat'\nalias qsa='cf --qstatall'" >> ~/.bashrc
```

> _Note:_ These tools don't work with LSF, as I don't have a LSF testing
> server to work on. Please get in touch if you can help.

### --qdel
Sometimes you may be running multiple pipelines and want to stop just one.
It can be a pain to find the job numbers to do this manually, so instead
you can use Cluster Flow to kill these jobs. When running `cf --qstat`,
ID values are printed for each pipeline. For example:
```
$ qs

======================================================================
 Cluster Flow Pipeline: fastq_bowtie
 Submitted:             17 hours, 1 minutes, 46 seconds ago
 Working Directory:     /path/to/working/dir
 Cluster Flow ID:       fastq_bowtie_1468357637
 Submitted Jobs:        29
 Running Jobs:          1
 Queued Jobs:           2 (dependencies)
 Completed Jobs:        26 (89%)
======================================================================
```

You can then use this _Cluster Flow ID_ to kill all jobs within that pipeline:

```
cf --qdel fastq_bowtie_1468357637
```

### --add_genome
Run the Cluster Flow interactive wizard to add new genomes.

### --setup
Run the interactive setup wizard to create a configuration file for
Cluster Flow.

### --version
Display the currently installed version of Cluster Flow.

### --check_updates
Check online for any available Cluster Flow updates.

### --help
Show a help message describing the different command line flags available.
