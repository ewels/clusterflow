---
title: General Usage
layout: toc
---

# Using Cluster Flow for the first time
This part of the documentation is intended for end users. If you're looking for how to set up Cluster Flow, see the [installation instructions]({{site.url}}/installation/)
## Environment Modules
Many compute clusters use [environment modules](http://modules.sourceforge.net/) to manage the user's `PATH`. This is entirely optional - if you're happy adding Cluster Flow to your `PATH` another way, just skip this section.

If you're not sure what all of that meant, try running `module avail clusterflow` on the command line. You should see the available module. If not, skip this section.

If you are using environment modules and Cluster Flow is set up to work with them, you'll need to load the Cluster Flow module before you can start using it. To do so, run the following command:

    module load clusterflow 

If you want to make your life easier when using environment modules, you can auto-load Cluster Flow each time you log in. To do this, run the following command:

    sed -i "$ a\module load clusterflow" ~/.bashrc

## Personal Config Files

### Config Wizard

Before you start using Cluster Flow, it’s a good idea to set yourself up with a personalised Cluster Flow config file. Cluster Flow has an interactive mode which will guide you through the required values and create a config file for you. To run this, run the following command:

	cf --make_config

### Genome Wizard

Most pipelines need a reference genome. Running `--add_genome` gives an interactive wizard which will lead you through the process of adding new genome paths. See the [Installation Instructions]({{site.baseurl}}/installation/#adding_genome_paths) for detailed instructions.

# Typical usage

To run cluster flow, use the `cf` command.

Most pipelines will need a genome name to be specified. IDs such as `NCBIM37` or `GRCh37` are used as keys to identify genome paths for aligning reads. These are specified by adding the flag `--genome` followed by an ID when you run Cluster Flow.

After flags are specified, the name of the module or pipeline to be run should be declared. If there is a pipeline with the same name as a module it will take preference.

Finally, write the input filenames. This can be done with linux wildcard expansion (`*`). For example, to run the `sra_bismark` pipeline, the command would be:

	cf --genome NCBIM37 sra_bismark *.sra

To read more about the command line options, see the [Command Line Reference]({{site.url}}/cl_reference/).

# Paired end / single end files
If using Cluster Flow with FastQ files, it will try to guess whether the files are paired end or single end. This is done by sorting the filenames, stripping `_[1-4]` and then comparing each file name to the next. If a pair is found to be identical, they are processed as paired end files.

When input files aren’t FastQ files, each Cluster Flow module will try to determine paird end or single end data by the same method. This behaviour is particularly useful when processing SRA files, as a single input file can split into multiple paired end files. Downstream modules can be run blindly and will behave correctly according to the output `.fastq` files produced by `fastq_dump`.
This behaviour can be overridden by specifying `--paired` or `--single` on the command line.

# Downloading files
When downloading files from an FTP server, parallelisation can be a problem. Cluster Flow has a download function built into its core `cf` package to deal with downloads. Cluster Flow sets off each download as a cluster job, with queue IDs set so that they process in series. The pipeline jobs also wait on the download IDs, so as soon as a download is finished it begins processing. In this way, processing does not need to wait for all downloads to finish, yet downloads are able to run in series and not overwhelm the server or internet connection.

To use this feature, simply submit download URLs instead of input filenames. Cluster Flow will recognise anything starting with `http`, `https` or `ftp` as a URL and set it off with the `download` module.

If using the `--file-list` parameter you can specify a filename to save each download. Add this on the same line as the download URL, separated by a tab character. This is particularly useful when downloading arbitrarily named SRA files. Labrador is able to generate these download files, offering a simple way to get download URLs and sensibly named files.

See below for an example download `--file-list` file:

	ftp://ftp-trace.ncbi.nlm.nih.gov/sra/SRR944695.sra Input_4OHT_rep5.sra
	ftp://ftp-trace.ncbi.nlm.nih.gov/sra/SRR944694.sra Input_4OHT_rep4.sra
	ftp://ftp-trace.ncbi.nlm.nih.gov/sra/SRR944693.sra Input_4OHT_rep3.sra
	ftp://ftp-trace.ncbi.nlm.nih.gov/sra/SRR944692.sra Input_4OHT_rep2.sra
	ftp://ftp-trace.ncbi.nlm.nih.gov/sra/SRR944691.sra Input_4OHT_rep1.sra
	
# Avoiding cluster overload
If using Cluster Flow with a large number of files it can be easy to swamp the available resources on the cluster and annoy any other users trying to get things done. Cluster Flow has several built in features to try to avoid this.

Firstly, CF can limit the number of parallel runs by using either `--split_files` or `--max_runs`. By default, Cluster Flow runs with a maximum of 12 parallel runs per pipeline. This default can be modified in the Cluster Flow config file or at run time.

In addition to limiting the number of parallel jobs it runs, Cluster Flow can try to limit the memory usage and number of cores each module uses. It calculates the maximum number of jobs that can be theoretically running at the same time when a pipeline is initiated. The available cores and memory are divided by this number and these ideal cores and memory per module are presented these numbers to each module. The modules can compare these values to their minimum and maximum requirements and request appropriate resources.

The end result of this process should be that jobs are run with maximum resources and speed when there are not too many files to process. Jobs will be run with minimum resources when processing many files so as not to overwhelm the available compute resources.

This behaviour is configured using the `@total_cores` and `@total_mem` configuration options, which can be overridden with the `--cores` and `--mem` command line arguments.