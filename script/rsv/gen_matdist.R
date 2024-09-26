#!/usr/bin/env Rscript
library(seqinr)
library(stringr)
library(ggplot2)
library(reshape2)
library(optparse)


option_list = list(
  make_option(c("--input"), type="character", default=NULL, 
              help="input aln file", metavar="character"),
  make_option(c("--matrix"), type="character", default=NULL, 
              help="output matrix pdf", metavar="character"),
  make_option(c("--table"), type="character", default=NULL, 
              help="output table tsv", metavar="character") 
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

allseq<-read.fasta(opt$input)

list_sample<-names(allseq)

#list_sample<-list_sample[1:4]
list_result<-list()
cmpt<-1
for(sample1 in list_sample){
  for(sample2 in list_sample){
    #sample1<-"22Pl116−022017870401HP_S470|WG"
    #sample2<-"22Pl116−022017870401HP_S470|WG"
    
    aln1<-toupper(allseq[[sample1]])
    aln2<-toupper(allseq[[sample2]])
    
    nb_diff<-0
    list_diff<-c()
    k<-1
    while(k<=length(aln1)){
      if((aln1[k]!="N" && aln2[k]!="N" && aln1[k]!=aln2[k])){
        newk<-k
        if(aln1[newk]=="-"){
          while(aln1[newk]=="-" & newk<=length(aln1)) newk<-newk+1
          list_diff<-c(list_diff,paste0("1:DEL",k,"-",newk))
          nb_diff<-nb_diff+1
          k<-newk
        } else if(aln2[newk]=="-"){
          while(aln2[newk]=="-" & newk<=length(aln2)) newk<-newk+1
          list_diff<-c(list_diff,paste0("2:DEL",k,"-",newk))
          nb_diff<-nb_diff+1
          k<-newk
        } else {
          list_diff<-c(list_diff,paste0("pos:",k,"/1:",aln1[k],"/2:",aln2[k]))
          nb_diff<-nb_diff+1
          k<-k+1
        }
      } else {
        k<-k+1
        newk<-k
        }
    }
    vec_result<-c(sample1,sample2,as.character(nb_diff),paste0(list_diff,collapse=";"))
    list_result[[cmpt]]<-vec_result
    cmpt<-cmpt+1
  }
}

table_merge<-as.data.frame(do.call("rbind", list_result))
colnames(table_merge)<-c("comp1","comp2","count","description")
table_merge$count<-as.numeric(as.character(table_merge$count))

write.table(table_merge,opt$table,sep="\t",row.names = F,quote=F)


matrix<-acast(table_merge,comp1~comp2,value.var="count")
matrix[lower.tri(matrix)]<-NA
table<-melt(matrix, na.rm = TRUE)
colnames(table)<-c("comp1","comp2","count")

final_plot<-ggplot(data = table, aes(comp1, comp2, fill = count))+
  xlab("")+ylab("")+
  geom_tile(color = "white",show.legend = FALSE)+
  geom_text(data=table,aes(comp1, comp2, label = count), color = "black", size = 6)+
  theme(axis.text.x = element_text(vjust = 1, 
                                   size = 6, hjust = 1,angle = 45),
        axis.text.y = element_text(vjust = 0.5, 
                                   size = 6, hjust = 1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 8))+
  scale_fill_gradient(low="blue",high = "red")

#final_plot

pdf(opt$matrix)
final_plot
dev.off()
