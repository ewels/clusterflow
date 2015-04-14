#!/usr/bin/env Rscript

#########################################
### CLUSTER FLOW COMMAND LINE ARGUMENTS
#########################################
args <- commandArgs(TRUE)
cf <- list()
cf[['params']] <- list()

for (i in 1:length(args) ){
    j <- i + 1
    if(substring(args[i], 0, 2) == '--'){
        # Next arg is also a flag or doesn't exist, this must be a setting
        if((exists(args[j]) && substring(args[j], 0, 2) == '--') || is.na(args[j])){
            cf[[substring(args[i],3)]] = TRUE
            next
        }

        # Next case not a flag, must be value
        else {

            # Special case - params
            if(args[i] == '--params'){
                vals <- unlist(strsplit(args[j], '='))
                cf[['params']][[vals[1]]] <- vals[2]
            }

            # Everything else
            else {
                cf[[substring(args[i],3)]] <- args[j]
            }

        }
    }
}

# Print help
if(!is.null(cf[['help']])){
    cat("------------\nmethylKit Cluster Flow Module\n------------\n")
    cat("This module takes sorted SAM files from bismark and calculates\n")
    cat("cross-sample correlations and plots some nice figures.\n")
    q(save="no")
}

# Requirements
if(!is.null(cf[['requirements']])){
    cat("cores:1\nmemory:10G\nmodules:\ntime:4:00:00\n")
    q(save="no")
}



######################
### SET UP
######################

# Running for real - check that we have everything
if(is.null(cf[['run_fn']]) || is.null(cf[['job_id']]) || is.null(cf[['prev_job_id']])){
    cat("\n###CF Error - Missing critical command line parameter (run_fn, job_id or prev_job_id)\n\n")
    print(args)
    q(save="no")
}

# Read the run file
file.list <- list()
runfile <- file(cf[['run_fn']], open="r")
rf <- readLines(runfile)
comment_block <- FALSE
for (i in 1:length(rf)){

    # Strip whitespace
    rf[i] <- gsub('^\\s+', '', rf[i])
    rf[i] <- gsub('\\s+$', '', rf[i])

    # Comments
    if(substring(rf[i],0,2) == '/*'){
        comment_block <- TRUE
        next
    }
    if(substring(rf[i],0,2) == '*/'){
        comment_block <- FALSE
        next
    }

    if(length(grep('^[@#>]', rf[i])) == 0 && comment_block == FALSE){
        sections <- unlist(strsplit(rf[i], "\t"))
        if(length(sections) > 0 && sections[1] == cf[['prev_job_id']]){
            file.list <- c(file.list, sections[2])
        }
    }
}
close(runfile)

if(length(file.list) == 0){
    cat("\n###CF Error - No files from previous job found!\n\n")
    q(save="no")
}

# Define the sample names
sample.list = file.list

# Make a fake treatment list - we don't know how the samples relate
treatments = as.list(rep(1, length(sample.list)))

######################
### INSTALLATION
######################

# Check if we have the methylKit module already
if("methylKit" %in% rownames(installed.packages()) == FALSE) {
  # Not found - install everything...

  # dependencies
  install.packages( c("data.table","devtools"))
  source("http://bioconductor.org/biocLite.R")
  biocLite(c("GenomicRanges","IRanges"))

  # install the development version from github
  library(devtools)
  install_github("al2na/methylKit", build_vignettes=FALSE)
}

# Load methylKit library
library(methylKit)


######################
### LOAD DATA
######################

# Load sorted deduplicated bismark BAM files to a methylRawList object
objs=read.bismark(location=file.list, sample.id=sample.list, assembly="GRCh37",
                  save.folder=NULL, save.context=NULL, read.context="CpG",
                  nolap=TRUE, mincov=10, minqual=20, treatment=treatments)
