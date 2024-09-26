#!/usr/bin/env Rscript
#v0.0.2
library(gtable)
library(grid)
library(gridExtra)
library(optparse)
library(cowplot)

option_list = list(
  make_option(c("--output"), type="character", default=NULL, 
              help="output file.", metavar="character"),
  make_option(c("--plot_folder"), type="character", default=NULL, 
              help="path to folder to be printed", metavar="character"),
  make_option(c("--chunk"), type="integer", default=8, 
              help="number of plot on one page (default=8)", metavar="character") 
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

lst_plot<-list.files(opt$plot_folder,pattern="*.RDS",full.names=TRUE)
lst_plot<-lapply(lst_plot,readRDS)
lst_plot<-lapply(lst_plot,as_grob)

pdf(opt$output,width=24,height=29,bg="white",colormodel="cmyk",paper="A4")
#nb : x plots are printed on each page, can be changed in chunk option
pc<-opt$chunk
i<-1
while (i <= length(lst_plot)){
  if (i+pc < length(lst_plot)){ #check if i+chunk reaches the end of the list
    lst_plot_sub<-lst_plot[i:(i+pc-1)]
    hts=rep(1,length(lst_plot_sub))
    grid.arrange(grobs=lst_plot_sub,ncol=1,nrow=opt$chunk,heights=c(hts))
  } else { #print remaining plot if i+chunk greater than list length
    print("last remaining")
    lst_plot_sub<-lst_plot[i:length(lst_plot)]
    hts=rep(1,length(lst_plot_sub))
    grid.arrange(grobs=lst_plot_sub,ncol=1,nrow=length(lst_plot_sub),heights=c(hts))
  }
  i<-i+pc
}
dev.off()