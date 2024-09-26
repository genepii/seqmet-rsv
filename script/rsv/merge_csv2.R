#!/usr/bin/env Rscript
#v0.0.1
library(optparse)
library(data.table)

option_list = list(
  make_option(c("--input"), type="character", default=NULL, 
              help="file listing input", metavar="character"),
  make_option(c("--output"), type="character", default=NULL, 
              help="output filename", metavar="character"),
  make_option(c("--separator"), type="character", default=",", 
              help="separator, default ','", metavar="character") 
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

files_list<-list.files(path=opt$input, pattern=".tsv|.csv", full.names=TRUE)

result <- rbindlist(lapply(files_list, function(x) fread(x, colClasses = 'character', sep = opt$separator, data.table = FALSE)), use.names=TRUE, fill=TRUE)
fwrite(result, opt$output)
