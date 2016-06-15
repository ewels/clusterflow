---
title: Usage
layout: toc
---

# General Usage
## Listing what's available
Once Cluster Flow is [up and running](installation), you can list the pipelines,
modules and reference genomes which are available using the following commands:

    cf --pipelines            # List available pipelines
    cf --modules              # List modules
    cf --genomes              # List reference genomes

## Getting help
To get instructions for how to use Cluster Flow on the command line, use:

    cf --help

You can also use this command to find out more information about
pipelines and modules:

    cf --help <module-name>
    cf --help <pipeline-name>

## Starting a run
In its most basic form, analyses are run as follows:

    cf <pipeline> <files>

Single modules can also be specified instead of a pipeline, _e.g._:

    cf <module> *.bam

Most pipelines and modules will need a reference genome, specified
using `--genome`:

    cf --genome GRCh37 sra_bowtie *.sra

The name following `--genome` is the name assigned when adding the reference
genome to Cluster Flow. This can be seen when listing genomes with `cf --genomes`.

To read more about the command line options, see the [reference page](reference).

## Filename checking
When launching Cluster Flow, a number of filename checks are performed. If given
FastQ files look like paired-end files, it launches in paired-end mode (this can
be overridden with `--single`).

If `@merge_regex` is configured in the configuration file, input files may be
merged before processing.

If a mixture of file types or paired end / single end FastQ files are found,
Cluster Flow will show an error and exit. This step can be skipped by using the
`--no_fn_check` parameter.

## Downloading files
As well as supplying Cluster Flow with input files, you can give URLs. This will
cause Cluster Flow to add the `cf_download` module to the start of your pipeline
to download the data.

Cluster Flow will recognise anything starting with `http`, `https` or `ftp` as
a URL. Downloads are processes in series to avoid overwhelming the internet connection.

If using the `--file-list` parameter you can also specify a filename for each download.
Add this on the same line as the download URL, separated by a tab character. This is
particularly useful when downloading arbitrarily named SRA files and is compatible
with the [Labrador Dataset Manager](https://github.com/ewels/labrador).

## Avoiding cluster overload
Cluster Flow has a number of features built in to avoid swamping your cluster with jobs.

Firstly, limits the number of parallel runs created. Defaults are set in the config file with
`@split_files` and `@max_runs`. These make Cluster Flow submit a default number of files per run
(usually one), but limits it to a maximum number of parallel runs (usually twelve).

Cluster Flow also try to intelligently limit the memory usage and number of cores each module uses.
The config options `@total_cores` and `@total_mem` specify the maximum resources to be used
by each Cluster Flow pipeline. These are split up amongst the max simultaneous jobs and
presented to each module. The modules can then request resources, making use of optional
parallelisation where available.
