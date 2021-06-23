library("readr")
library(stringr)
#Extraction of data
args <- commandArgs(trailingOnly=TRUE)

IDs <- read_tsv(args[1], skip=0,col_names=FALSE)

IDs.table <- as.data.frame(matrix(ncol=length(t(IDs)),nrow=0))
names(IDs.table)=as.vector(t(IDs))

write_tsv(IDs.table,"Genes_ID-base.txt",col_names=TRUE)
