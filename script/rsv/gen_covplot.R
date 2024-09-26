#!/usr/bin/env Rscript
#v0.0.2
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(gtable)
library(grid)
library(gridExtra)
library(optparse)
library(egg)

option_list = list(
  make_option(c("--output"), type="character", default=NULL, 
              help="output file in .RDS format", metavar="character"),
  make_option(c("--bga_file"), type="character", default=NULL, 
              help="bga bed formated file .", metavar="character"),
  make_option(c("--sample_name"), type="character", default=NULL, 
              help="name of the sample.", metavar="character"), 
  make_option(c("--ref_name"),type="character", default=NULL,
              help="name of the ref", metavar="character")
);
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

### Function Section ###

#function that create a gtable from a list of ggplot
get_ggplot<-function(lst_ggplot,lst_size_seg,name){
  panel=c()
  ht=c()
  whole_name=c(name)
  #calculate width and heigth for each panel
  for(size in lst_size_seg){
    percent<-size/max(unlist(lst_size_seg))
    panel=c(panel,percent)
    ht=c(ht,0.5)
  }
  whole_name<-c(whole_name,rep.int("", length(lst_size_seg)-1))
  #this function create a gtable with panels containing each segment
  final_plot<-egg::ggarrange(plots=lst_ggplot,ncol=length(lst_size_seg),widths=panel,labels=whole_name,label.args=list(gp = grid::gpar(fontsize = 10),hjust=-0.5))
  return(final_plot)
}

#function to create a ggplot, base on the number of segment
sub_bga<-function(segment,idx,ref){
bga_sub<-subset(bga, bga$segment==x)
if(ref=="MN908947"){
    if(idx>1){ #case where no legend is required (except for x axis)
      plot<- ggplot(bga_sub) + 
        geom_rect(aes(xmin=start,xmax=end,ymin=0,ymax=depth),color="blue") +
        xlab(segment)+
        ylim(0,log10(100000))+
        geom_hline(yintercept=log10(10),color="red")+
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank(), 
              axis.title.y=element_blank())
    }
    else{ #case where legend is required
      plot<- ggplot(bga_sub) + 
        geom_rect(aes(xmin=start,xmax=end,ymin=0,ymax=depth),color="blue") +
        xlab(segment)+
        ylab("log10 depth")+
        ylim(0,log10(100000))+
        geom_hline(yintercept=log10(10),color="red")+
        geom_vline(xintercept = 21563, linetype="dotted",color="black",linewidth=0.8)+
        geom_vline(xintercept = 25384,linetype="dotted",color="black",linewidth=0.8)+
        theme()
      }
    }
    else {
      if(idx>1){ #case where no legend is required (except for x axis)
        plot<- ggplot(bga_sub) + 
          geom_rect(aes(xmin=start,xmax=end,ymin=0,ymax=depth),color="blue") +
          xlab(segment)+
          ylim(0,log10(100000))+
          geom_hline(yintercept=log10(10),color="red")+
          theme(axis.text.y=element_blank(),
                axis.ticks.y=element_blank(), 
                axis.title.y=element_blank())
      }
      else{ #case where legend is required
        plot<- ggplot(bga_sub) + 
          geom_rect(aes(xmin=start,xmax=end,ymin=0,ymax=depth),color="blue") +
          xlab(segment)+
          ylab("log10 depth")+
          ylim(0,log10(100000))+
          geom_hline(yintercept=log10(10),color="red")+
          theme()
      }
    }
  return(plot)
}

#open file
bga_raw<-read.table(opt$bga_file,header=FALSE,sep="\t")
names(bga_raw)<-c("chromosome","start","end","depth")

bga_raw$start<-as.numeric(bga_raw$start)
bga_raw$end<-as.numeric(bga_raw$end)
bga_raw$depth<-as.numeric(bga_raw$depth)
#convert depth in log scale 
bga_raw$depth<-log10(bga_raw$depth+1)

#set max_depth
bga_raw$depth<-ifelse(bga_raw$depth>=log10(100000),log10(100000),bga_raw$depth)

#create new col for segment contained in the "chromosome" name
bga<-separate(data = bga_raw, col = chromosome, into = c("left", "right"), sep = "\\|")
names(bga)<-c("ref","segment","start","end","depth")
#remove NA for REF
bga$segment[is.na(bga$segment)]<-"REF"
#get uniq segment name
seg<-unique(bga$segment)

lst_size_seg=list()
lst_ggplot=list()

#create x variables each contaning only their segment
i=1
for (x in seg){
  #store segment's plot into list
  lst_ggplot[[i]]<-sub_bga(x,i,opt$ref_name)
  #get length of each segment
  bga_seg<-subset(bga, bga$segment==x)
  lst_size_seg[[i]]<-max(bga_seg$end)
  i=i+1
}

### we have : lst of ggplot for every segment, length for every segment, sequence name
### we want : create a gtable that regroup every segment
title <-paste0(opt$ref_name," - ",opt$sample_name)
gtable<-get_ggplot(lst_ggplot,lst_size_seg,title)

#export the gtable in a RDS format
saveRDS(gtable,file=opt$output)