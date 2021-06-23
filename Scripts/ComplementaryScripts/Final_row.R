library('readr')
library('stringr')
#Extraction of data
args=commandArgs(trailingOnly=TRUE)

file <- paste("Genes",args[1],".txt",sep='')

df<-as.data.frame(read_tsv(file,col_names=TRUE))

All <- nrow(df)+1
HL <- nrow(df)+2
LL <- nrow(df)+3

df[All:LL,1]<-0 

df[is.na(df)]<-0

#Calculation of the absolute number of CNVs (both duplication and deletions)
  for (i in colnames(df)){
    df[All,i]<-sum(abs(df[[i]]))
    df[HL,i]<-sum(abs(df[1:5,i]))
    df[LL,i]<-sum(abs(df[6:10,i]))
  }

name_csv=paste("Final_table",args[1],".csv",sep='')
name_tsv=paste("Final_table",args[1],".txt",sep='')
write_csv(df,name_csv,col_names=TRUE)
write_tsv(df,name_tsv,col_names=TRUE)

#rownames(df) <- c("HL0421","HL0438","HL0626","HL0672","HL0677","LL0409","LL0703","LL0720","LL0733","LL1010","Total_CNVs","HL_CNVs","LL_CNVs")


