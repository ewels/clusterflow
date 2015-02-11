Cluster Flow
============

Cluster Flow is a pipelining tool to automate and standardise bioinformatics analyses on high-performance cluster environments. It is designed to be easy to use, quick to set up and flexible to configure.

## Cluster Flow Website - [http://ewels.github.io/clusterflow/](http://ewels.github.io/clusterflow/)

There's a new website which for the Cluster Flow documentation which has loads of helpful information and examples. You can see it here: [http://ewels.github.io/clusterflow/](http://ewels.github.io/clusterflow/)

If you're anxious to just get your hands on the code, check out the [releases page](https://github.com/ewels/clusterflow/releases)

Licence
-------
Cluster Flow is released with a GPL v3 licence. Cluster Flow is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. For more information, see the licence that comes bundled with Cluster Flow.

Change Log
----------
#### v0.4 devel
* **Warning: Break of backwards compatability**
	* The way that genome references are handled has been rewritten. They are now agnostic
      to the type of reference, making it far easier to whatever type of reference you need.
      Additionally, the wizard to add genome paths has been written and is now largely automated,
      making it super fast to add new genomes.
	* A consequence of this change is any `genomes.config` files written before v0.3 of
      Cluster Flow will no longer work. Thankfully the fix is easy! Replace `@bowtie_path`
      with `@reference  bowtie`. `@gtf_path` changes to `@reference gtf` and so on.
      `@genome_path` changes to `@reference  fasta`.
	* If you have any custom pipelines these will also need to be updated. `@require_bowtie`
      changes to `@require_reference bowtie`and so on. See updated example module files
      for examples on how to update custom modules.
	* Apologies for any inconvenience that this change incurs. Feel free to [get in touch](https://github.com/ewels)
      if you have any problems.
* New Stuff
	* Support for STAR RNA-seq aligner (thanks to [@stu2](https://github.com/stu2))
	* Modules are given more information via the run file to help
       decide the amount of memory and cores they bid for (eg. number of files, reference)
	* All perl scripts now have `env perl` in shebang to increase portability
	* Now using GRIDEngine `h_vmem` memory option instead of `vf`
        * Gives a hard memory limit instead of a request limit at job submission time
        * Thanks to  [@stu2](https://github.com/stu2) and [@s-andrews](https://github.com/s-andrews)
	* Support for explicit GRIDEngine queue nomination on the command line
	* New `--verbose` and `--print_versions` flags
        * Prints version information where possible. Supported by bowtie 1 & 2, STAR, samtools_sort, trim_galore, tophat and 
            htseq-count modules. This information is sent to the report file.
        * Prints the commands submitted to the cluster to `STDOUT` when launched.
    * New `--runfile_prefix` option to help avoid potential filename clashes.
    * Added compatability with GRIDEngine `~/.sge_request files` (by ignoring them). 
      Thanks to [@s-andrews](https://github.com/s-andrews)
    * New tophat module called `tophat` which introduces a workaround for buggy MAPQ 
        reporting by tophat whilst keeping unique alignments. Thanks to [@FelixKrueger](https://github.com/FelixKrueger).
        * The previous tophat module is still available if you're not interested in MAPQ scores and 
            would like slightly faster processing. This is now called `tophat_broken_MAPQ.cfmod`.
    * Pipeline completion e-mails are now written to disk as well (HTML and plain text)
    * New log file containing the job submission commands as well as the output received from the cluster at submission (usually numeric job identifiers)
    * Removed the `--qstatcols` command line option and added the `@colourful` config option to replace it
        * The config wizard is also updated to add this to your personal config
* Bugs Squashed
	* Fixed output filename problem in tophat with output cleaning
	* Fixed bugs causing minimum memory allocation regardless of availability
	* Fixed bug causing Bowtie2 to fail if Bowtie1 index absent
    * Cleaned up some unrecognised output that always made it into the log file

#### [v0.3](https://github.com/ewels/clusterflow/releases/tag/v0.3) - 2014-07-11
* New Stuff
	* Awesome new HTML report e-mails
		* Much more readable HTML report e-mails which look super-snazzy (see [example](http://ewels.github.io/clusterflow/example_report_good.html))
		* Any errors are highlighted making them quick to identify (see [example](http://ewels.github.io/clusterflow/example_report_bad.html))
		* Custom strings set in the config can flagged as [highlights](http://ewels.github.io/clusterflow/example_report_highlights.html) or as [warnings](http://ewels.github.io/clusterflow/example_report_warnings.html)
		* Designed to work on desktop and mobile phone screens
	* Cluster Flow now re-orders the log file so that output from different modules doesn't overlap
		* Made each module prepend its stdout and sterr with a CF module flag
		* Made the `cf_run_finished` module parse the above flag and print out module by module
	* Rewrote how the environment module loading works - now much more robust
		* Uses `modulecmd` to prepare perl syntax commands
		* Moved code into the `load_environment_modules` function in Helpers.pm
	* Added environment module aliases
		* This allows you to load specific environment module versions or use different names to those specified within CF modules
		* eg. Replace `fastqc` with `FastQC/0.11.2`
* Updates
	* Made the `bismark_methXtract` module create genome-wide coverage reports if GTF path is available
	* Added the `-q` parameter to the FastQC module to make the log files cleaner
	* Removed the now uneccesary `bismark_tidy` module and renamed `bismark_messy` to `bismark_report`
* Bugs Squashed
	* Fixed dependency bug introduced in v0.2 which was making all downloads fire simultaneously
	* Fixed issue where modules using the CF::Constants Perl Module couldn't load the central config file
	* Fixed typo in environment module loading in `bismark_align` module
	* Reordered loading of the environment modules in `trim_galore` so that FastQC is loaded first, fixing dependency issues
	
#### [v0.2](https://github.com/ewels/clusterflow/releases/tag/v0.2) - 2014-05-29
* New Stuff
	* Now compatable with SLURM
	* Customise batch job commands in the config (see the [docs](http://ewels.github.io/clusterflow/installation/#making_cluster_flow_work_with_your_environment)
	* Created new GitHub pages website to hold documentation: http://ewels.github.io/clusterflow
* Updates
	* Ported repository to github: https://github.com/ewels/clusterflow
	* Wrote new readme for github
* Bugs Squashed
	* Custom modules in `~/clusterflow/modules/` weren't being found
	* General code clean-ups all over the place 

#### [v0.1](https://github.com/ewels/clusterflow/releases/tag/v0.1) - 2014-04-25
* The first public release of Cluster Flow, although it's been in use at the Babraham Institute for around 6 months. It's been in heavy development throughout that time and is now approaching a state of being relatively stable.


Credits
-------
Cluster Flow was written by [Phil Ewels](http://phil.ewels.co.uk) whilst working in the [Babraham Bioinformatics](http://www.bioinformatics.babraham.ac.uk/) group in Cambridge, UK. He now maintains it whilst working at the [Science for Life Laboratory](http://www.scilifelab.se/) in Stockholm, Sweden.

Cluster Flow has also had contributions from [@stu2](https://github.com/stu2), [@s-andrews](https://github.com/s-andrews) and [@FelixKrueger](https://github.com/FelixKrueger), amongst others.
