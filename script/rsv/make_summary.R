#!/usr/bin/env Rscript
#v0.1.0
library(optparse)
library(tidyr)
library(seqinr)
library(stringr)


option_list = list(
  make_option(c("--output"), type="character", default=NULL, 
              help="output file files.", metavar="character"),
  make_option(c("--refname"), type="character", default=NULL, 
              help="name of the reference considered.", metavar="character"),
  make_option(c("--coverage"), type="character", default=NULL, 
              help="coverage data from all_cov.cov file.", metavar="character"),
  make_option(c("--fastqstat-table"), type="character", default=NULL, 
              help="count data from qc_fastqstat.tsv file.", metavar="character"),
  make_option(c("--bamstat-table"), type="character", default=NULL, 
              help="count data from all_bamstat.tsv file.", metavar="character"),
  make_option(c("--count-table"), type="character", default=NULL, 
              help="count data from all_counts.tsv file.", metavar="character"),
  make_option(c("--coinf-table"), type="character", default=NULL, 
              help="coinf data from coinf tsv file.", metavar="character"),
  make_option(c("--cons"), type="character", default=NULL, 
              help="fasta file containing all generated consensus.", metavar="character"),
  make_option(c("--posc"), type="character", default=NULL, 
              help="covidseq positive controls data from posc.tsv file.", metavar="character"),
  make_option(c("--conta"), type="character", default=NULL, 
              help="contamination results from contamination_common_poolt.tsv file.", metavar="character"),
  make_option(c("--veriftable"), type="character", default=NULL, 
              help="premap_verif.tsv file.", metavar="character"),
  make_option(c("--runname"), type="character", default=NULL, 
              help="run name", metavar="character"),
  make_option(c("--reference"), type="character", default=NULL, 
              help="reference for denovo analysis", metavar="character"),
  make_option(c("--blast"), type="character", default=NULL, 
              help="blast results from denovo analysis", metavar="character")   
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

checkPolySeg<-function(reference_id){
  split<-str_split(reference_id,coll("|"))
  if(length(split[[1]])>1) result<-as.character(split[[1]][[2]]) else result<-"FAIL"
  return(result)
}

getReference<-function(reference_id){
  split<-str_split(reference_id,coll("|"))
  result<-as.character(split[[1]][[1]])
  return(result)
}

print("#consensus")
## FASTA
if(!(is.null(opt$cons))){
  alnfa <- read.fasta(file = as.character(opt$cons))
  seqnames<-names(alnfa)
  lengthseq = unlist(lapply(alnfa, function(l) length(l)))
  length_N_tot = unlist(lapply(alnfa, function(l) length(grep("n",l))))
  #Sequence trimming for the internal N check
  alnfa_trim <- alnfa
  for (sequence in names(alnfa_trim)){
    consensus<-as.vector(unlist(alnfa_trim[[sequence]]))
    trimcons<-consensus
    alnfa_trim[[sequence]]<-as.list(trimcons)
  }
  Nseqtrim = unlist(lapply(alnfa_trim, function(l) length(grep("n",l))))
  N_check_results<-as.data.frame(cbind(seqnames,length_N_tot,lengthseq,Nseqtrim))
  rownames(N_check_results)<-NULL
  N_check_results$sample_id<-sapply(N_check_results$seqnames,getReference)
  N_check_results$seqtype<-sapply(N_check_results$seqnames,checkPolySeg)
  #if(lengthseq>0){
  if(length(readLines(opt$cons))>=2 && lengthseq>0){
    N_check_results$consensus_perccoverage<-100-(as.numeric(as.character(N_check_results$length_N_tot))/as.numeric(as.character(N_check_results$lengthseq))*100)
  }else{
    N_check_results$consensus_perccoverage<-"NA"
  }
  nb_seqtype<-length(unique(as.character(N_check_results[which(N_check_results$seqtype!="FAIL"),]$seqtype)))
  #check if polyseg
  N_check_results<-N_check_results[,-c(which(colnames(N_check_results)=="length_N_tot"),
                                       which(colnames(N_check_results)=="Nseqtrim"),
                                       which(colnames(N_check_results)=="lengthseq"))]
  if(nb_seqtype>1){
    N_check_results<-pivot_wider(N_check_results,id_cols = "sample_id",
                                 names_from = "seqtype",
                                 values_from = c("consensus_perccoverage"))
    
    if(length(grep("FAIL",colnames(N_check_results)))>0) N_check_results<-N_check_results[,-grep("FAIL",colnames(N_check_results))]
    colnames(N_check_results)[2:length(colnames(N_check_results))]<-paste0("consensus_perccoverage_",colnames(N_check_results)[2:length(colnames(N_check_results))])
  }else{
    N_check_results<-N_check_results[,c("sample_id","consensus_perccoverage")]
  }
}

#reference compute contig size
## FASTA
if(!(is.null(opt$reference))){
  alnfa <- read.fasta(file = as.character(opt$reference))
  reference_id<-names(alnfa)
  reference_id<-unique(sapply(reference_id,function(x) str_split(x,coll("|"))[[1]][[1]]))
  if(length(reference_id)>1){
    stop(reference_id)
  }
  reference_length = sum(unlist(lapply(alnfa, function(l) length(l))))
  contigsize_table<-as.data.frame(cbind(reference_id,reference_length))
}

print("#coverage")
if(!(is.null(opt$coverage))){
  cov_table<-read.table(paste0(as.character(opt$repres),as.character(opt$coverage)),header=F) 
  colnames(cov_table)<-c("reference_id","seqtype","bam_meandepth","sample_id")
  cov_table$polyseg<-sapply(cov_table$seqtype,checkPolySeg)
  nb_seqtype<-length(unique(as.character(cov_table[which(cov_table$seqtype!="FAIL"),]$seqtype)))
  #check if polyseg
  if(nb_seqtype>1){
    mean_cov<-pivot_wider(cov_table,id_cols = c("sample_id","reference_id"),names_from = "polyseg",values_from = "bam_meandepth")
    if(length(grep("FAIL",colnames(mean_cov)))>0) mean_cov<-mean_cov[,-grep("FAIL",colnames(mean_cov))]
    colnames(mean_cov)[3:length(colnames(mean_cov))]<-paste0("bam_meandepth_",colnames(mean_cov)[3:length(colnames(mean_cov))])
  }else{
    mean_cov<-cov_table[,c("sample_id","reference_id","bam_meandepth")]
  }
}

print("#blast")
if(!(is.null(opt$blast))){
  if(file.info(opt$blast)$size>300){
    tableblast<-read.delim(opt$blast,header=T,colClasses = "character",sep="\t")
    refname<-opt$refname
    tableblast$reference_id<-refname
    tableblast$QUERY_SEQID<-sapply(tableblast$QUERY_SEQID,function(x) str_split(x,coll("|"))[[1]][[1]])
    #keep best results using indentity percentage and bitscore
    maxbit<-aggregate(as.numeric(tableblast$BITSCORE),by=list(tableblast$QUERY_SEQID),max)
    colnames(maxbit)<-c("QUERY_SEQID","MAX_BITSCORE")
    tableblast<-merge(tableblast,maxbit,by="QUERY_SEQID")
    tableblast<-subset(tableblast,as.numeric(BITSCORE)==MAX_BITSCORE)
    maxperc<-aggregate(as.numeric(tableblast$PERC_IDENT),by=list(tableblast$QUERY_SEQID),max)
    colnames(maxperc)<-c("QUERY_SEQID","MAX_PERC_IDENT")
    tableblast<-merge(tableblast,maxperc,by="QUERY_SEQID")
    tableblast<-subset(tableblast,as.numeric(PERC_IDENT)==MAX_PERC_IDENT)
    tableblast<-tableblast[with(tableblast, order(QUERY_SEQID)),]
    #get only one result among the good ones
    tableblast$nbres<-""
    for(currsample in unique(tableblast$QUERY_SEQID)) tableblast[tableblast$QUERY_SEQID==currsample,]$nbres<-seq(1,nrow(tableblast[tableblast$QUERY_SEQID==currsample,]))
    tableblast<-subset(tableblast,nbres==1)
    tableblast<-tableblast[c("QUERY_SEQID","reference_id","SUBJECT_TITLE","BITSCORE","PERC_IDENT")]
    colnames(tableblast)<-c("sample_id","reference_id","blast_subjectid","blast_bitscore","blast_percident")
  } else {
    opt$blast<-NULL
  }
}

print("#verif")
if(!(is.null(opt$veriftable))){
  verif_table<-read.table(opt$veriftable,header=F,colClasses = "character")
  alnfa <- read.fasta(file = as.character(opt$reference))
  reference_id<-names(alnfa)
  segmentnb<-length(reference_id)
  colnames(verif_table)<-as.character(verif_table[1,])
  verif_table<-verif_table[-1,]
  verif_table<-pivot_longer(
  data=verif_table,
  cols=all_of(colnames(verif_table)[2:length(colnames(verif_table))]),
  names_to = "sample_id",
  values_to = "mapping")
  verif_table$mapping<-as.numeric(verif_table$mapping)
  verif_table$reference_id<-sapply(verif_table$SAMPLEID,function(x) str_split(x,coll("|"))[[1]][[1]])
  verif_table$segment<-sapply(verif_table$SAMPLEID,function(x) str_split(x,coll("|"))[[1]][[2]])
  verif_table$mapsup90<-ifelse(verif_table$mapping>=0.9,1,0)

  #get max mapping for each segment
  maxsegverif_table<-aggregate(verif_table$mapping,by=list(verif_table$sample_id,verif_table$segment),max)
  colnames(maxsegverif_table)<-c("sample_id","segment","maxmapping")

  verif_table<-merge(verif_table,maxsegverif_table,by=c("sample_id","segment"))
  verif_table$ismax<-ifelse(verif_table$mapping==verif_table$maxmapping,1,0)
  maxrefverif_table<-aggregate(verif_table$ismax,by=list(verif_table$sample_id,verif_table$reference_id),sum)
  colnames(maxrefverif_table)<-c("sample_id","reference_id","nbsegmax")
  maxrefverif_table$count<-ifelse(maxrefverif_table$nbsegmax==segmentnb,1,0)

  uniqrefverif_table<-aggregate(maxrefverif_table$count,by=list(maxrefverif_table$sample_id),sum)
  colnames(uniqrefverif_table)<-c("sample_id","count")
  uniqrefverif_table$asuniqmaxref<-ifelse(uniqrefverif_table$count==1,"OK","FAILED")
  uniqrefverif_table<-uniqrefverif_table[,c("sample_id","asuniqmaxref")]

  #count segments above 0.9 for all reference
  segverif_table<-aggregate(verif_table$mapsup90,by=list(verif_table$sample_id,verif_table$segment),sum)
  colnames(segverif_table)<-c("sample_id","segment","count")
  segverif_table<-aggregate(segverif_table$count,by=list(segverif_table$sample_id),sum)
  colnames(segverif_table)<-c("sample_id","count")
  segverif_table$hassegmultiverif<-ifelse(segverif_table$count>segmentnb,"WARNING","OK")
  segverif_table<-segverif_table[,c("sample_id","hassegmultiverif")]

  finalverif_table<-merge(uniqrefverif_table,segverif_table,by="sample_id")

  #final statements
  print(finalverif_table)
  finalverif_table$bam_verif<-""
  finalverif_table$bam_verif<-ifelse(finalverif_table$asuniqmaxref=="OK" & finalverif_table$hassegmultiverif=="OK","OK",finalverif_table$bam_verif)
  finalverif_table$bam_verif<-ifelse(finalverif_table$asuniqmaxref=="FAILED","FAILED",finalverif_table$bam_verif)
  finalverif_table$bam_verif<-ifelse(finalverif_table$asuniqmaxref=="OK" & finalverif_table$hassegmultiverif=="WARNING","WARNING",finalverif_table$bam_verif)
  finalverif_table<-finalverif_table[,c("sample_id","bam_verif")]
}

print("#fastqstat")
if(!(is.null(opt$fastqstat))){
  fastqstat_table<-read.table(as.character(opt$fastqstat),header=F, sep="\t")
  colnames(fastqstat_table)<-c("reference_id","type","count","sample_id")
  fastqstat_table<-pivot_wider(fastqstat_table,id_cols = c("reference_id", "sample_id"),names_from = "type",values_from = "count")
  fastqstat_table$reference_id<-gsub("RAW", opt$refname, fastqstat_table$reference_id)
}

print("#bamstat")
if(!(is.null(opt$bamstat))){
  bamstat_table<-read.table(as.character(opt$bamstat),header=F, sep="\t")
  colnames(bamstat_table)<-c("reference_id","type","count","sample_id")
  bamstat_table<-bamstat_table[order(bamstat_table$count, decreasing=TRUE),]
  bamstat_table<-pivot_wider(bamstat_table,id_cols = c("reference_id","sample_id"),names_from = "type",values_from = "count")
}

print("#count")
if(!(is.null(opt$count))){
  count_table<-read.table(as.character(opt$count),header=F, sep="\t")
  colnames(count_table)<-c("reference_id","type","count","sample_id")
  count_table<-pivot_wider(count_table,id_cols = c("reference_id","sample_id"),names_from = "type",values_from = "count")
}

print("#coinf")
if(!(is.null(opt$coinf))){
  coinf_table<-read.table(as.character(opt$coinf),header=F, sep="\t")
  colnames(coinf_table)<-c("reference_id","type","value","sample_id")
  coinf_table<-pivot_wider(coinf_table,id_cols = c("reference_id","sample_id"),names_from = "type",values_from = "value")
}

print("#posc")
if(!(is.null(opt$posc))){
  posc_header<-read.table(as.character(opt$posc), nrows=1, header=F, stringsAsFactors=FALSE)
  posc_table<-read.table(as.character(opt$posc), skip=1, header=F)
  posc <- data.frame(sample_id = character(0), qc_seqcontrol = character(0), stringsAsFactors=FALSE)
  for (i in 2:length(posc_header)) {posc[nrow(posc)+1,]<-c(posc_header[i], any(posc_table[,i] >= 0.05))}
  posc$qc_seqcontrol<-gsub("FALSE", "FAILED", posc$qc_seqcontrol)
  posc$qc_seqcontrol<-gsub("TRUE", "OK", posc$qc_seqcontrol)
}

print("#conta")
if(!(is.null(opt$conta))){
  conta_header<-read.table(as.character(opt$conta), nrows=1, header=F, stringsAsFactors=FALSE)
  conta_table<-read.table(as.character(opt$conta), skip=1, header=F)
  conta <- data.frame(sample_id = character(0), vcf_dpcount = character(0), stringsAsFactors=FALSE)
  #for (i in 2:length(conta_header)) {conta[nrow(conta)+1,]<-c(conta_header[i], as.character(conta_table[,i][which.max(lapply(conta_table[,i], function(x) str_split(x, "/", n=Inf, simplify = TRUE)[,1]))]))}
  for (i in 2:length(conta_header)) {conta[nrow(conta)+1,]<-c(conta_header[i], as.character(lapply(conta_table[,i][which.max(lapply(conta_table[,i], function(x) str_split(x, "/", n=Inf, simplify = TRUE)[,1]))], function(x) str_split(x, "/", n=Inf, simplify = TRUE)[,1])))}
}

print("#merge")
alldata<-fastqstat_table
print(head(alldata))
if(!(is.null(opt$bamstat))) alldata<-merge(alldata,bamstat_table,by = c("sample_id","reference_id"), all.x=TRUE)
print(head(alldata))
if(!(is.null(mean_cov))) alldata<-merge(alldata,mean_cov,by = c("sample_id","reference_id"), all.x=TRUE)
print(head(alldata))
#refname not found in N_check_results, TODO
if(!(is.null(N_check_results))) alldata<-merge(alldata,N_check_results,by = "sample_id", all.x=TRUE)
print(head(alldata))
if(!(is.null(opt$reference))) alldata<-merge(alldata,contigsize_table,by = "reference_id", all.x=TRUE)
print(head(alldata))
if(!(is.null(opt$conta))) alldata<-merge(alldata,conta,by = "sample_id", all.x=TRUE) else {
  alldata$vcf_dpcount<-NA
}
if(!(is.null(opt$posc))) alldata<-merge(alldata,posc,by = "sample_id",all.x=TRUE)
print(head(alldata))
if(!(is.null(opt$veriftable))) alldata<-merge(alldata,finalverif_table,by="sample_id",all.x=TRUE)
print(head(alldata))
if(!(is.null(opt$blast))) alldata<-merge(alldata,tableblast,by = c("sample_id","reference_id"), all.x=TRUE) else {
  alldata$blast_subjectid<-NA
  alldata$blast_bitscore<-NA
  alldata$blast_percident<-NA
}
print(head(alldata))
if(!(is.null(opt$count))) alldata<-merge(alldata,count_table,by = c("sample_id","reference_id"), all.x=TRUE)
print(head(alldata))
if(!(is.null(opt$coinf))) alldata<-merge(alldata,coinf_table,by = c("sample_id","reference_id"), all.x=TRUE)
print(head(alldata))
if(!(is.null(opt$runname))) alldata$run_id<-as.character(opt$runname)

#write.csv2(alldata,as.character(opt$output),quote=F,row.names = F)
write.table(alldata,as.character(opt$output),quote=F,row.names = F,sep="\t")
