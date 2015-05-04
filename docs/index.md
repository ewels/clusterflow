---
title: Cluster Flow
layout: default
---

# Cluster Flow

<p class="lead">A command-line pipeline tool which uses common cluster managers to run bioinformatics analysis pipelines.</p>

Benefits of using Cluster Flow:

* Routine analyses are very quick to run, for example: `cf --genome GRCh37 fastq_bowtie *fq.gz`
* Pipelines use identical parameters, standardising analysis and making results more reproducable
* Integrated parallelisation tools help prevent your cluster becoming overloaded
* All commands and output is logged in files for future reference
* Intuitive commands and a comprehensive manual make Cluster Flow easy to use
* Job monitoring tools and E-mail notifications allow you to keep track of your jobs

How Cluster Flow differs from other pipeline tools:

* Very lightweight and flexible
* Pipelines and configurations can easily be generated on a project-specific basis if required
* New modules and pipelines are super easy to write (see video tutorial)

Cluster Flow currently supports GRIDEngine (SGE), LSF and SLURM, plus it should be fairly easy to port to others.

<div class="demo_gif inactive">
<img src="_files/CF_mini.gif" title="mini CF demo">
</div>

## Tutorial Videos
There are two tutorial videos up on YouTube, see below:

* [Usage / Installation Tutorial](http://youtu.be/b2g_zQiz9ys) - How to configure and run Cluster Flow
* [Advanced Tutorial](http://youtu.be/aBHOcsA2M6w) - How to write your own pipelines and modules

Please see the sidebar for access to the full Cluster Flow documentation. If you're looking for a quick overview of what Cluster Flow is and how it works, see the [Introduction](/introduction/).

## Contributing to Cluster Flow
If you write a module or pipeline which could be of use to others, or modify the core Cluster Flow code in a helpful way, it would be great to merge those changes back into the core Cluster Flow project.

The easiest way to do this is to [fork the Cluster Flow repository](https://help.github.com/articles/fork-a-repo), make your changes, committing them and pushing them as you go. When you've finished, submit a [pull request](https://help.github.com/articles/using-pull-requests) and the new code can be merged into the central project.

## Change Log

<h4 class="version"><a href="https://github.com/ewels/clusterflow/releases/tag/v0.3"><span class="label label-success version">Version 0.3</span></a> <small>2014-07-11</small></h4>

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


<h4 class="version"><a href="https://github.com/ewels/clusterflow/releases/tag/v0.2"><span class="label label-success version">Version 0.2</span></a> <small>2014-05-29</small></h4>

* New Stuff
	* Now compatable with SLURM
	* Customise batch job commands in the config (see [docs]({{site.baseurl}}/installation/#making_cluster_flow_work_with_your_environment) for more info)
	* Created new GitHub pages website to hold documentation: [http://ewels.github.io/clusterflow/](http://ewels.github.io/clusterflow/)
* Updates
	* Ported repository to github: [https://github.com/ewels/clusterflow](https://github.com/ewels/clusterflow)
	* Wrote new readme for github
* Bugs Squashed
	* Custom modules in `~/clusterflow/modules/` weren't being found
	* General code clean-ups all over the place

<h4 class="version"><a href="https://github.com/ewels/clusterflow/releases/tag/v0.1"><span class="label label-success version">Version 0.1</span></a> <small>2014-04-25</small></h4>

* The first public release of Cluster Flow, although it's been in use at the Babraham Institute for around 6 months. It's been in heavy development throughout that time and is now approaching a state of being relatively stable.
