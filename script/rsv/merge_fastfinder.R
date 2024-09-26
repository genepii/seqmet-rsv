library(tidyr)
library(stringr)
library(optparse)
#v0.0.4

option_list = list(
  make_option(c("--output"), type="character", default=NULL, 
              help="output file name.", metavar="character"),
  make_option(c("--val_input"), type="character", default=NULL, 
              help="*_validation.tsv separated by a \",\" ", metavar="character") 
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

output_file<-as.character(opt$output)
#output_file<-"test_fastfinder_merged.csv"

#read input files
if(!(is.null(opt$val_input))){
  input_string<-as.character(opt$val_input)
  #setwd("C:/Users/regueha/Desktop/")
  #input_string<-"000000_H3N2_validation.tsv,000000_pH1N1_validation.tsv"
  input_list<-as.character(str_split(input_string,",")[[1]])
  files_list<-list()
  for(files in input_list) {
    table<-read.delim(files,sep="\t",header=T,colClasses = "character")
    table<-table[,c("sample_id","val_glims","summary_reference_id","val_result","val_commentary")]
    files_list[[files]]<-table
    }
}

#regroup table
table_merge<-do.call("rbind", files_list)
rownames(table_merge)<-NULL

#keep only necessary columns
#table_merge<-table_merge[,c("val_glims","summary_reference_id","val_result","val_commentary")]
#write.csv(table_merge,"tablemerge.csv")

#table_merge<-unique(table_merge)

#get both "val_result" and "val_commentary" values for each reference in one single row per sample
table_pivot<-pivot_wider(table_merge,
                         id_cols=c("sample_id"),
                         names_from="summary_reference_id",
                         values_from =c("val_result","val_commentary") )

write.csv2(table_pivot,"test_pivot.csv")

#for each sample, if one unique val_result value (i.e ININT) for all reference.
#if not, remove the ININT values and keep the unique one instead.
table_pivot$GABILL<-""
GABILL_table<-table_pivot[,grep("val_result",colnames(table_pivot))]
#for each sample, if one unique val_commentary value (i.e DSC) for all reference.
#if not, remove the DSC values and keep the unique one instead.
table_pivot$COMGRAB<-""
COMGRAB_table<-table_pivot[,grep("val_commentary",colnames(table_pivot))]

for(i in 1:nrow(table_pivot)){
  #i<-1
  #val_result_vec<-as.character(table_pivot[i,grep("val_result",colnames(table_pivot))])
  val_result_vec<-unlist(as.list(GABILL_table[i,]))
  print(val_result_vec<-unlist(as.list(GABILL_table[i,])))
  #if not REPASSE / ININT
  if(length(subset(val_result_vec,!(val_result_vec %in% c("ININT","REPASSE SEQ FAILED","NA"))))!=0){ 
    #if("TRUE" %in% grepl("COINF-A",val_result_vec) && "TRUE" %in% grepl("COINF-B",val_result_vec)){
    #if("TRUE" %in% grepl("RSVA|RSVB",val_result_vec) && "TRUE" %in% grepl("RSVA|RSVB",val_result_vec)){
    #  table_pivot[i,which(colnames(table_pivot)=="GABILL")]<-"COINF VRS A/B"
    #  table_pivot[i,which(colnames(table_pivot)=="COMGRAB")]<-"COINF VRS A/B"
    #  next
    #}
    table_pivot[i,which(colnames(table_pivot)=="GABILL")]<-subset(val_result_vec,!(val_result_vec %in% c("ININT","REPASSE SEQ FAILED")))
    #get commentary for selected reference result
    reference<-colnames(table_pivot)[which(table_pivot[i,]==subset(val_result_vec,!(val_result_vec %in% c("ININT","REPASSE SEQ FAILED"))))]
    reference<-str_split(reference,"_")[[1]][[3]]
    table_pivot[i,which(colnames(table_pivot)=="COMGRAB")]<-val_commentary_vec_f<-table_pivot[i,paste0("val_commentary_",reference)]
    # if REPASSE > ININT
  } else if("REPASSE SEQ FAILED" %in% val_result_vec){
    table_pivot[i,which(colnames(table_pivot)=="GABILL")]<-"REPASSE SEQ FAILED"
    table_pivot[i,which(colnames(table_pivot)=="COMGRAB")]<-"ININT"
    # ININT
  }else {
    table_pivot[i,which(colnames(table_pivot)=="GABILL")]<-"ININT"
    table_pivot[i,which(colnames(table_pivot)=="COMGRAB")]<-"DSC"
  }
}

#format results to fastfinder format
table_pivot<-table_pivot[,c("sample_id","GABILL","COMGRAB")]

table_pivot$sample_id<-sapply(table_pivot$sample_id,function(x) str_split(x,"-")[[1]][[2]])
table_pivot$sample_id<-sapply(table_pivot$sample_id,function(x) str_split(x,"_")[[1]][[1]])

fastfinder_table<-pivot_longer(table_pivot,cols=c("GABILL","COMGRAB"),names_to = "AssayResult")
fastfinder_table$instrument<-"NB552333"
fastfinder_table$author<-"laurence.josset@chu-lyon.fr"
fastfinder_table<-fastfinder_table[,c("sample_id","instrument","author","AssayResult","value")]
colnames(fastfinder_table)<-c("Sample ID","Instrument ID(s)","Analysis authorized by","AssayResultTargetCode","AssayResult")

write.csv(fastfinder_table,output_file,row.names = F,quote=F)
