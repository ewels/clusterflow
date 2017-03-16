## Listing what's available
Once Cluster Flow is up and running, you can list available pipelines,
modules and reference genomes which are available using the following commands:
```bash
cf --pipelines            # List pipelines
cf --modules              # List modules
cf --genomes              # List reference genomes
```

## Getting help
To get instructions for how to use Cluster Flow on the command line, use:
```bash
cf --help
```

You can also use this command to find out more information about
pipelines and modules:
```bash
cf --help [module-name]
cf --help [pipeline-name]
```

## Starting a run
In its most basic form, analyses are run as follows:
```bash
cf [pipeline] [files]
```

Single modules can also be specified instead of a pipeline:
```bash
cf [module] *.bam
```

Most pipelines and modules will need a reference genome, specified
using `--genome`:
```bash
cf --genome GRCh37 sra_bowtie *.sra
```

The ID following `--genome` is the ID assigned when adding the reference
genome to Cluster Flow. This can be seen when listing genomes with `cf --genomes`.

## Module Parameters
The default execution of different tools can be modified by using module
_parameters_. These can be set within pipeline scripts or on the command line.
Specifying `--param [example]` will apply the `[example]` parameter to every
module in the pipeline.

Different module support different parameters. Some are flags, some are key pairs.
To find out more, see the Modules documentation.

Typical things you can do are to set adapter trimming preferences with TrimGalore!:
```bash
cf --genome GRCh37 --param clip_r1=6 --param min_readlength=15 sra_bowtie *sra
```

or run Bismark in PBAT mode:
```cf
cf --genome GRCm38 --param pbat fastq_bismark *.fastq.gz
```

When setting in pipeline scripts, simply add the paramters after the module
names (tab-delimited). For example, this is the `trim_bowtie_miRNA` pipeline:
```
#trim_galore adapter=ATGGAATTCTCG
	#bowtie mirna
```
This sets a custom adapter for trimming and tells the bowtie module to use
the `mirna` parameter.


## Filename checking
When launching Cluster Flow, a number of filename checks are performed. If input
files are FastQ and the filenames look like paired-end files, it launches in
paired-end mode (this can be overridden with `--single`).

If a mixture of file types or paired end / single end FastQ files are found,
Cluster Flow will show an error and exit. This step can be skipped by using the
`--no_fn_check` parameter.

If `@merge_regex` is configured in the configuration file, matching input files
will be merged before processing.

## Downloading files
As well as supplying Cluster Flow with input files, you can give URLs. This will
cause Cluster Flow to add the `cf_download` module to the start of your pipeline
to download the data.

Cluster Flow will recognise anything starting with `http`, `https` or `ftp` as
a URL. Downloads are processes in series to avoid overwhelming the internet connection.

If using the `--file-list` parameter you can also specify a filename for each download.
This should be added after the download URL, separated by a tab character. This is
particularly useful when downloading arbitrarily named SRA files and is compatible
with the [Labrador Dataset Manager](https://github.com/ewels/labrador). See also
the stand-alone [SRA-Explorer](http://ewels.github.io/sra-explorer/) tool if you
don't have Labrador installed.

## Avoiding cluster overload
Cluster Flow has a number of features built in to avoid swamping your cluster
with jobs.

Firstly, Cluster Flow limits the number of parallel runs created. Defaults are
set in the config file with `@split_files` (default 1) and `@max_runs` (default 12).
`@split_files` defines the minimum number of files per run, `@max_runs` defines
the maximum number of parallel runs and adds more files per run if needed.

Cluster Flow also try to intelligently limit the memory usage and number of
cores each module uses. The config options `@total_cores` and `@total_mem`
specify the maximum resources to be used by each Cluster Flow pipeline. These
are split up amongst the max simultaneous jobs and presented to each module.
The modules can then request resources, making use of optional parallelisation
where available.
