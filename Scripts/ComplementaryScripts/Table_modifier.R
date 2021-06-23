library("readr")
library(stringr)
#Extraction of data
args=commandArgs(trailingOnly=TRUE)

name <- str_replace(args[1],"_ID-base",args[3])
if (args[4] == 0){
	IDs <- as.data.frame(read_tsv(args[1],skip=0,col_names=TRUE))
} else {
	IDs <- as.data.frame(read_tsv(name,skip=0,col_names=TRUE))
}

Indv <- as.data.frame(read_tsv(args[2],skip=0,col_names=TRUE))
Name <- unlist(strsplit(args[2],"/"))[2] 
new_R <- nrow(IDs)+1
IDs[new_R,1]<-0 
rownames(IDs[new_R,]) <- Name

rownames(Indv) <- Indv[["Gene"]]

for (i in Indv[["Gene"]]){
	if (Indv[i,"Type"] == "Dup"){
		IDs[new_R,i]=1
        }else if (Indv[i,"Type"] == "Del"){
		IDs[new_R,i]=-1
	}
}
IDs[is.na(IDs)]<-0
write_tsv(IDs,name,col_names=TRUE)