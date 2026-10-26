#!/usr/bin/env Rscript
file_arg <- grep('^--file=', commandArgs(), value=TRUE)
if (length(file_arg) != 1L) stop('Invoke with Rscript /path/to/run.R')
setwd(dirname(normalizePath(sub('^--file=', '', file_arg))))
source('code/pipeline.R')
pipeline_main(commandArgs(trailingOnly=TRUE))
