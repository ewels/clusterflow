---
title: Troubleshooting
layout: toc
---

# Bugs and Errors
If you come across a strange looking error message or find a bug, please do let us know. You submit new issues here: [https://github.com/ewels/clusterflow/issues](https://github.com/ewels/clusterflow/issues "Browse and submit new issues about Cluster Flow on github")

## Feature Requests
If you'd like Cluster Flow to do something it doesn't, log a request! The issue tracker system mentioned above can be used for enhancement requests too.

## Going old-school
If you can't face all of this fancy-pants github stuff, feel free to drop the author an e-mail at <a href="mailto:phil.ewels@scilifelab.se">phil.ewels@scilifelab.se</a>


# Frequently Asked Questions

## General

### Permission Errors
A number of errors can be caused by scripts not having executable file privileges. You can see the file permissions with `ls -l`, you should see something like this:

	$ ls -l clusterflow/modules/
	total 608
	-rwxrwxr-x 1 phil phil 6770 May 20 16:40 bismark_align.cfmod
	-rwxrwxr-x 1 phil phil 4291 May 16 12:54 bismark_deduplicate.cfmod
	-rwxrwxr-x 1 phil phil 2748 May 16 12:54 bismark_messy.cfmod
	-rwxrwxr-x 1 phil phil 5652 May 16 12:54 bismark_methXtract.cfmod
	-rwxrwxr-x 1 phil phil 3553 May 16 12:54 bismark_tidy.cfmod
	-rwxrwxr-x 1 phil phil 7119 May 16 12:54 bowtie1.cfmod

This example is for the modules directory (all modules should have executable privileges for all), the same applies to the main `cf` file.

### DOS carriage returns
If you've edited any files, you may get problems due to windows-based editors putting DOS-style `\r` carriage returns in.

Most linux environments come with a package called `dos2unix` which will clean these up. eg:
	
	dos2unix *


## Errors from environment modules

### ERROR:105: Unable to locate a modulefile for 'clusterflow'
This error probably means that Cluster Flow isn't installed in your environment module system, and you're trying to run `module load clusterflow`

You can skip this step if you have another way of accessing the `cf` file, or see the [Installation Instructions]({{site.baseurl}}/installation/#environment_modules) for details about how to set Cluster Flow up with environment modules.

## Errors from job submission

### Unable to run job: job rejected: the requested parallel environment "orte" does not exist.
This message means that your GRIDEngine setup doesn't have the default `orte` environment set up.  If you have different environments set up you can list them with:

	qconf -spl

You can get the details of any environment with:

	qconf -sp [name]

If you find one which assigns slots to a single node (`allocation_rule` should be `$fill_up`) then you can just do a find &amp; replace for `orte` to the name of your local environment and that should make things work again.

_(Answered by Simon Andrews)_

### Unable to run job: job rejected (other reasons)
There may be other differences in the job submission requests that cause them to fail. If you see errors such as this, you can use the `@custom_job_submit_command` configuration variable to customise the way that jobs are requested.

For more information, see the [Installation Instructions]({{site.baseurl}}/installation/#making_cluster_flow_work_with_your_environment)