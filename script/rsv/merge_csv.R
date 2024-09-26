#!/usr/bin/env Rscript
#v0.0.1
library(optparse)
library(data.table)
library(openxlsx)
library(stringi)

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

subsString <- function(string, corresTable, oldColumn = 1, newColumn = 2){
	corres <- as.character(corresTable[,newColumn][which(corresTable[,oldColumn] %in% string)])
	
	if(length(corres)==0){
		return(string)
	}
	else{
		return(corres)
	}
}

catColumns <- function(dataframe, mode = '', name, columnsList){
	colList<-unlist(stri_split_fixed(columnsList, ","))
	colList<-colList[colList %in% colnames(dataframe)]
	oldCol<-dataframe[, colList]
	oldCol<-lapply(oldCol, function(x) stri_replace_all_regex(stri_replace_all_regex(x, "([-]+)$", ""), "([.]+)$", ""))
	newCol<-do.call(paste, c(oldCol, sep = "|"))
	
	newCol<-lapply(newCol, function(x) stri_replace_all_fixed(x, ",", "."))
	newCol<-sapply(newCol, function(x) reduceCat(x, mode = mode))
	
	newDf<-data.frame(newCol)
	colnames(newDf)<-name
	
	return(cbind(dataframe, newDf))
}

reduceCat <- function(string, mode = ''){
	if(mode=='pcr'){
		if('NEG' %in% stri_split_regex(string, "[|]")[[1]] || 'ININT' %in% stri_split_regex(string, "[|]")[[1]] || '{<NONREAL}' %in% stri_split_regex(string, "[|]")[[1]] || '{<ININTERP}' %in% stri_split_regex(string, "[|]")[[1]]){
			cat<-c('99.9')
		}
		else{
			values<-lapply(stri_split_fixed(string, "|")[[1]], function(x) stri_replace_all_regex(stri_replace_all_fixed(x, "SARS-CoV2", "NCOV"), "[^0-9.]", ""))
			values<-values[sapply(values, stri_length) > 0]
			cat<-do.call(paste, c(values, sep = ","))
		}
		if(length(cat)==0){cat<-NA}
	}
	else{
		values<-lapply(stri_split_fixed(string, "|")[[1]], function(x) stri_replace_all_regex(stri_replace_all_regex(x, "^-+", ""), "-+$", ""))
		values<-values[sapply(values, stri_length) > 0]
		cat<-do.call(paste, c(values, sep = ","))
		if(length(cat)==0){cat<-NA}
	}
	
	return(cat)
}

glimsXlsx <- function(path, firstRow = 12, columnsCat = data.frame(), namesCol = 1, namesColSimplify = FALSE, rowNamesSubs = data.frame()){
	df <- readWorkbook(xlsxFile = path, sheet = 1, startRow = firstRow, colNames = TRUE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE, skipEmptyCols = TRUE, rows = NULL, cols = NULL, check.names = FALSE, sep.names = "_", namedRegion = NULL, na.strings = "", fillMergedCells = FALSE)
	
	if(namesColSimplify){
		rownames(df)<-lapply(df[,namesCol], function(x) stri_replace_all_fixed(stri_split_regex(x, "[()]")[[1]][1], " ", ""))
	}
	else{
		rownames(df)<-df[,namesCol]
	}
	
	if(length(rowNamesSubs) != 0){
		rownames(df)<-lapply(rownames(df), function(x) subsString(x, namesSubs))
	}
	
	if(length(columnsCat) != 0){
		for (i in 1:length(columnsCat[,1])) {
			df<-catColumns(df, as.vector(columnsCat[i,][1]), as.vector(columnsCat[i,][2]), as.vector(columnsCat[i,][3]))
		}
	}
	
	return(df)
}

catFile<-read.table(as.character(opt$catfile), header=F, sep="\t", encoding='UTF-8')
catFile<-sapply(catFile, function(x) stri_replace_all_fixed(x, " ", "_"))

df1 <- readWorkbook(xlsxFile = "N_VD_CONF_RESPI_HOL_IAI_cribinfo_190601-200601.xlsx", sheet = 1, startRow = 12, colNames = TRUE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE, skipEmptyCols = TRUE, rows = NULL, cols = NULL, check.names = FALSE, sep.names = "_", namedRegion = NULL, na.strings = "", fillMergedCells = FALSE)
rownames(df1)<-df1[,1]

result <- rbindlist(lapply(readLines(opt$input), function(x) fread(x, colClasses = 'character', sep = opt$separator, data.table = FALSE)), use.names=TRUE, fill=TRUE)
fwrite(result, opt$output)
