## Overview
Modules are the heart of Cluster Flow. Each module is a wrapper around
a single bioinformatics tool. Each module has three modes of operation:

1. Specifying the required resources for the job
2. Running the bioinformatics tool
3. Printing a help message about the module

Modules are executed using system commands, so can be written in any
language. However, most existing modules are written in Perl.

Module filenames must be in the format `<module_name>.cfmod.<extension>`,
eg. `mymod.cfmod.pl`. They can be stored in the following locations
(chosen in this order of preference):

* Current working directory
* `~/.clusterflow/modules/`
* `<installation_dir>/modules/`

### Example module
An example module comes bundled with Cluster Flow, containing some highly
commented pseudocode which you can modify for your own uses. You can see
it in your `modules` directory:
[`example_module.pl`](https://github.com/ewels/clusterflow/blob/master/modules/example_module.pl)

### Existing Perl Scripts
If you have an existing script or tool, it's tempting to try to convert
it into a Cluster Flow module. However, I recommend instead keeping it
as a standalone script and creating a Cluster Flow module to launch this
instead. In our experience, this is much easier. It also has the
advantage that your script can still be run outside Cluster Flow.

## Specifying resources
At the top of every Cluster Flow module is a hash that defines the resources
needed by the tool. It looks something like this:

```perl
my %requirements = (
    'cores' 	=> $cores,
    'memory' 	=> $mem,
    'modules' 	=> $modules,
    'references'=> $refs,
    'time' 		=> $time
);
```

Each of these variables can be specified as a string, an array specifying
a range of appropriate values, or a subroutine to calculate a value based
on information specific to the run.

### Cores
The number of required cores can be specified either as a string or an array.
If your tool always uses a fixed number of cpus (for example, 1 if it's not
multi-threaded), just specify that number in quotes (`'cores' => '1'`).

If your tool can be sped up by using multiple cpus, you can specify a minimum
and maximum number in an array (`'cores' => ['3','8']`). Cluster Flow will
then allocate a number within that range according to how many jobs are being
created in parallel. This way, jobs will run as fast as possible for a handful
of files, but not overwhelm the cluster if many are being run at once.

### Memory
Memory works just like cores, above - either specify a string or an array with
a minimum and maximum amount. Numbers with no suffix will be interpreted
as bytes, then you can use `K`, `M` and `G` suffixes to specify kilobytes, megabytes
and gigabytes (`'memory' => '8G'`).

In some cases it can be useful to use a subroutine to dynamically calculate
the required memory. For example, you could inspect the filesize of a fasta
genome reference to determine the required memory:
```perl
'memory' 	=> sub {
    my $cf = $_[0];
    if (defined($cf->{'refs'}{'fasta'}) && -e $cf->{'refs'}{'fasta'}) {
        # Multiple the reference filesize (in bytes) by 1.2
        my $mem_usage = int(1.2 * -s $cf->{'refs'}{'fasta'});
        return CF::Helpers::bytes_to_human_readable($mem_usage);
    } else {
        # Sensible default
        return '8G';
    }
},
```

### Modules
A string or array of strings describing environment modules that should be loaded.
Try to keep this as generic as possible. People can specify specific versions or
naming in personal config files using `@environment_module_alias`.

### References
Genome reference and annotation is labelled with a field to describe it's type.
If a reference is required, you should specify its type here. This prevents Cluster
Flow from being launched if the reference genome is not specified.

For example, the bowtie2 module specifies `'references'=> 'bowtie2'`; the
featureCounts module specifies `'references' => 'gtf'`.

### Time
Some HPC clusters require a time limit to be specified when launching jobs. Here
you should predict approximately how long your module should run.

Some modules will always take a fixed amount of time to run, in which case this
can be specified as a string. For ten minutes, specify `'time' => '10'`.

The execution time for most modules will depend on how many input files they
are processing. **Modules often run with multiple sets of input files.**
To cope with this, supply a subroutine to this variable which can flexibly
request an amount of time according to how many input files will be processed.

The helper function `minutes_to_timestamp` is useful here - it takes a number
of minutes and returns a properly formatted timestamp (see below for more
information about helper functions).

If a module typically takes three hours to run, it could request it as follows:
```perl
'time' => sub {
    my $cf = $_[0];
    my $num_files = $cf->{'num_starting_merged_aligned_files'};
    return CF::Helpers::minutes_to_timestamp ($num_files * 3 * 60);
}
```
The `$cf` variable is a hash containing information about the job. See below
for a description of the keys available.

Remember to be conservative - high time requests can delay queue priority,
but low time requests will result in job failure.

### Help Text
Cluster Flow can request help text from a module if called with
`cf --help <module_name>`. You should write some text describing what the
module does, including any parameters or customisation available.

```perl
my $helptext = "".("-"x15)."\n My awesome module\n".("-"x15)."\n
This module is brilliant and worked first time because the author
read all of the Cluster Flow documentation! What a hero!\n\n";
```

## Module launch
Once the requirements hash and help text are written, we call a core
helper function called `module_start`. If the module is being called to
request resource requirements or help, the function will exit. If it is
being executed in a cluster job, it will return as hash with useful
information such as the input filenames.

Requirements should be passed to the function as a reference:
```perl
my %cf = CF::Helpers::module_start(\%requirements, $helptext);
```
The returned hash contains the following keys:
_(**NB**: Not all of these are available in request subroutines)_

```perl
%cf = {
    refs = '<hash>',                # Reference annotation for the specified genome. Keys are the reference type, values are the path to the annotation.
    prev_job_files = '<array>',     # File names resulting from preceding job.
    starting_files = '<array>',     # File names for the initial files that this thread of the pipeline was started with.
    files = '<hash>',               # Hash of arrays with all files from this pipeline thread. Keys are the module job IDs, values are arrays of output files.
    cores = '<int>',                # The number of cores allocated to the module.
    memory = '<str>',               # The amount of memory allocated to the module.
    params = '<hash>',              # A hash of key: value pairs. Value is `True` if only a flag.
    config = '<hash>',              # Hash containing arbitrary key: value configuration pairs from the run file. Always contains hash with key `notifications`.
    num_starting_files = '<int>',   # Number of files that this thread of the pipeline started with.
    num_starting_merged_files = '<int>', # Number of files that this thread of the pipeline started with, after merging if matched merge regex.
    num_starting_merged_aligned_files = '<int>',  # Guess at number of files after alignment, based on whether pipeline is running in paired end mode or not.
    pipeline_id = '<str>',          # The unique Cluster Flow ID of this pipeline. Useful for generating filenames.
    pipeline_name = '<str>',        # The name of the pipeline that was launched.
    pipeline_started = '<int>',     # A unix timestamp of when the pipeline was started.
    job_id = '<str>',               # The unique Cluster Flow ID for this job.
    prev_job_id = '<str>',          # The unique Cluster Flow ID for the previous job in the pipeline.
    run_fn = '<str>',               # The filename of the run file for this thread of the pipeline.
    run_fns = '<array>',            # All run file filenames for this pipeline (summary modules only).
    modname = '<str>',              # Name of this module
    mod_fn = '<str>',               # Filename of this module
}
```

### Checks
Although not necessary, most modules that use genome references do a sanity
check to make sure that they have what they need after this point. For
example, the STAR module checks that it has the required reference:
```perl
# Check that we have a genome defined
if(!defined($cf{'refs'}{'star'})){
	die "\n\n###CF Error: No star path found in run file $cf{run_fn} for job $cf{job_id}. Exiting..";
} else {
	warn "\nAligning against $cf{refs}{star}\n\n";
}
```

### Version logging
Again not necessary, but good practice - modules log the version of software
that they're about to run for future reference:
```perl
warn "---------- < module > version information ----------\n";
warn `MY_COMMAND --version`;
warn "\n------- End of < module > version information ------\n";
```

### Parameters
Modules are able to customise the way that they run depending on the presence
of custom parameters are run time. These are used for a range of reasons, such
as customising bowtie alignments for miRNA data, changing trimming settings
depending on library preparation type and many others. You can basically use
them however you like, though you'll find may modules doing this sort of thing:
```perl
my $extra_flag = (defined($cf{'params'}{'myflag'})) ? '--extra_flag' : '';
my $specific_var = '';
if(defined($cf{'params'}{'myvar'})){
    $specific_var = '--myvar '.$cf{'params'}{'myvar'};
}
# ..later..
$cmd = "mycommand --always $extra_flag $specific_var"
```

### Opening the run file
Each part of the pipeline has a `.run` file, used by the modules to track
the configuration options and output filenames as the pipeline progresses.

Your pipeline will need to open this run file in append mode so that it can
add the file names of any output that it creates.
```perl
open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";
```

## Command execution
### Looping through files
Once you have everything ready, you'll want to actually run your tool.
Remember that modules typically run with a collection of input files, so you
will need to loop through these and process them in sequence.

How you do this looping depends on what input your tool expects. If your
tool takes a single file and doesn't care whether it's paired end or single,
you can simply loop through all files from the previous job:
```perl
foreach my $file (@{$cf{'prev_job_files'}}){
    # process $file
}
```

Most preprocessing and alignment tools need either one single end FastQ
file or two paired end FastQ files. To handle this, you can use the
`is_paired_end` helper function to separate the input files into single end
and paired end:
```perl
my ($se_files, $pe_files) = CF::Helpers::is_paired_end(\%cf, @{$cf{'prev_job_files'}});
```
These files can then be looped over in separate loops:
```perl
# Go through each single end files and run Bowtie
if($se_files && scalar(@$se_files) > 0){
	foreach my $file (@$se_files){
        # process $file
    }
}
if($pe_files && scalar(@$pe_files) > 0){
	foreach my $files_ref (@$pe_files){
		my @files = @$files_ref;
        if(scalar(@files) == 2){
            # process $files[0] and $files[1]
        } else {
			warn "\n###CF Error! Bowtie paired end files had ".scalar(@files)." input files instead of 2\n";
		}
    }
}
```

### Building a command
Typically, Cluster Flow modules build a system command in a string. This is
then printed to stderr with the `###CFCMD` prefix. This is picked up by
Cluster Flow and added to the summary html report and e-mail.
```perl
my $command = "my_command -i $file -o $output_fn";
warn "\n###CFCMD $command\n\n";
```

### Running the command
Once build, the command should be executed using the perl `system` command.
This command returns the exit code once complete, which can be checked to
see whether the module has worked or not:
_(0 is success, which evaluates to false)_
```perl
if(!system ($command)){
    # command worked
} else {
    # Command returned a non-zero result, probably went wrong...
    warn "\n###CF Error! Example module (SE mode) failed for input file '$file': $? $!\n\n";
}
```

### Adding the output to the run file
If your command ran successfully, you should have created a new output
file. This should be added to the `.run` file along with the current
job id, so that it can be used by subsequent modules in the pipeline:
```perl
if(-e $output_fn){
    print RUN "$cf{job_id}\t$output_fn\n";
} else {
    warn "\n###CF Error! Example module output file $output_fn not found..\n";
}
```

## Job Completion
### Run File Output
A run file is created by Cluster Flow for each batch of files. It describes
variables to be used, the pipeline specified and the filenames used by each
module. The syntax of variables and pipeline is described in Pipeline syntax.

File names are described by a job identifier followed by a tab then a filename.
Each module is provided with its own job ID and the ID of the job that was run
previously. By using these identifiers, the module can read which input files
to use and write out the resulting filenames to the run file when complete.
Example run file syntax:
```
first_job_938 filename_1.txt
first_job_938 filename_2.txt
second_job_375 filename_1_processed.txt
second_job_375 filename_2_processed.txt
```

There can be any number of extra parameters, these are specific the module
and are specified in the pipeline configuration.

### E-mail report highlights
Any `STDOUT` or `STDERR` that your module produces will be written to the
Cluster Flow log file. At the end of each run and pipeline, an e-mail will
be sent to the submitter with details of the run results (if specified by
the config settings). Because the log file can be very long Cluster Flow pulls
out any lines starting with `###CF`. Typically, such a line should be printed
when a module finishes, with a concise summary of whether it worked or not.
Messages including the word `Error` will be highlighted and cause the final
e-mail to have warning colours. The configuration options `@log_highlight_string`
and `@log_warning_string` can customise this reporting.

Modules should print the command that they are going to run to STDERR so that
this is recorded in the log file. These are also sent in the e-mail
notification and should start with `###CFCMD`.

### Exit codes
It's likely that your cluster will continue to fire off the dependent jobs as
soon as the parent jobs finish, irrespective of their output. If a module
fails, the cleanest way to exit is with a success code, but without printing
any resulting output filename. The following modules will not find their
input filenames and so should immediately exit.


## Appendices
### Command line flags
Cluster Flow modules are expected to respond to the following command
line flags:

Flag                  | Step | Description
----------------------|------|---------------------------------------
`--requirements`      | 1    | Request the cluster resources needed by the module
`--run_fn` <str>      | 2    | Path to the Cluster Flow run file(s) for this pipeline
`--job_id` <str>      | 2    | Cluster job ID for this job
`--prev_job_id` <str> | 2    | Cluster job ID for the previous job
`--cores <int>`       | 2    | Number of cores allocated to the module
`--mem <str>`         | 2    | Amount of memory allocated to the module
`--param <str>`       | 2    | Extra parameters to be used
`--help`              | 3    | Print module help

The step number refers to whether the module is being executed:

1. By the core Cluster Flow script at pipeline launch
2. Within a cluster job, executing the tool
3. By the core Cluster Flow script, when `cf --help <modname>` is specified.


## Helper Functions
If your module is written in Perl, there are some common Cluster Flow
packages (perl modules) that you can use to provide some pre-written
functions.

There are currently three packages available to Cluster Flow modules.
`Helpers` contains subroutines of general use for most modules.
`Constants` and `Headnodehelpers` contain subroutines primarily for
use in the main `cf` script. You can include the Helpers package by
adding the following to the top of your module file:
```perl
use FindBin qw($Bin);
use lib "$FindBin::Bin/../source";
use CF::Helpers;
```

We use the package `FindBin` to add the binary directory to the path
(where `cf` is executing from).

Note that there is a Python version of the Helpers script which contains
many of the same functions and works in a comparable way.

### module_start
Handles the initiation of all modules. See above for a description of use.

### parse_runfile
Parses `.run` files. Called by `module_start()` and not usually run directly.

### load_environment_modules
Used to load environment modules into the PATH. Typically used by the
main `cf` script, though occasionally used elsewhere for special occasions.

### is_paired_end
This function takes an array of file names and returns an array of single
end files and an array of arrays of paired end files.

First, it checks the configuration set in the `.run` file.
If `@force_paired_end` is set, it sorts the files from the last job into
pairs and returns them. If `@force_single_end` is set it returns all
previous files as single end.

If neither variables are set, it sorts the files alphabetically, then
removes any occurance of `_[1-4]` from the filename and compares the list.
Identical pairs are returned as paired end.
```perl
my ($se_files, $pe_files) = CF::Helpers::is_paired_end(@$files);
foreach my $file (@$se_files){
    print "$file is single end\n";
}
foreach my $files_ref (@$pe_files){
    my @files = @$files_ref;
    print $files[0]." and ".$files[1]." are paired end.\n";
}
```

### is_bam_paired_end
Looks at BAM/SAM file headers and tries to determine whether it has been
generated using paired end input files or single end. The function reads
through the first 1000 reads of the file and counts how many `0x1` flags
it finds (denoting a paired read). If ``>= 800` of those first 1000 reads
are paired end, it returns `true`.
```perl
if(CF::Helpers::is_bam_paired_end($file)){
    ## do something with paired end BAM
} else {
    ## do something with single end BAM
}
```

### fastq_encoding_type
Scans a FastQ file and tries to determine the encoding type. Returns strings
`integer`, `solexa`, `phred33`, `phred64` or `0` if too few reads to safely
determine. This is done by observing the minimum and maximum quality scores.

For more details, see the
[Wikipedia page on FastQ encoding](http://en.wikipedia.org/wiki/FASTQ_format#Encoding)
```perl
($encoding) = CF::Helpers::fastq_encoding_type($file);
```

### fastq_min_length
Scans the first 100000 reads of a FastQ file and returns the longest read
length that it finds.
```perl
my $min_length = CF::Helpers::fastq_min_length($file);
```

### parse_seconds
Takes time in seconds as an input and returns a human readable string.
The optional second `$long` variable determines whether to use `h`/`m`/`s`
(`0`, false) or `hours`/`minutes`/`seconds` (`1`, true, the default).
```perl
my $time = CF::Helpers::parse_seconds($seconds, $long);
```

### timestamp_to_minutes / minutes_to_timestamp
Functions to convert between a SLURM / HPC style timestamp and minutes.
Attempts string parsing in the following order:

1. minutes
2. minutes:seconds
3. hours:minutes:seconds
4. days-hours
5. days-hours:minutes
6. days-hours:minutes:seconds

```perl
my $minutes = CF::Helpers::timestamp_to_minutes($timestamp);
my $timestamp = CF::Helpers::minutes_to_timestamp($minutes);
```

### human_readable_to_bytes / bytes_to_human_readable
Two functions which convert between human readable memory strings
(eg. `4G` or `100M`) and bytes.
```perl
my $bytes = CF::Helpers::human_readable_to_bytes('3G');
my $size = CF::Helpers::bytes_to_human_readable('7728742');
```

### mem_return_mbs
Takes a memory string and returns a number of megabytes, rounding
up to the nearest MB.

### allocate_cores
Takes the suggested number of cores to use, a minimum and maximum number
and returns a sensible result.
```perl
my $cores = CF::Helpers::allocate_cores($recommended, $min, $max);
```

### allocate_memory
Takes the suggested number of memory to use, a minimum and maximum amount
and returns a sensible result. Input can be human readable strings or bytes.
Returns a value in bytes.
```perl
my $mem = CF::Helpers::allocate_memory($recommended, $min, $max);
```

### cf_compare_version_numbers
Function to properly compare software version numbers. Correctly returns
that `v0.10` is greater than `v0.9`.

