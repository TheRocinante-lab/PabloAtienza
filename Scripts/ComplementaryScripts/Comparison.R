library("readr")
library(stringr)
#Extraction of data
args=commandArgs(trailingOnly=TRUE)

GenesData <- read_tsv(args[1], skip = 0, col_names = c("Chr","Start","End","Gene","Mean"))
if ("Relation" %in% colnames(GenesData)){
	GenesData$Relation <- NULL
	if ("Type" %in% colnames(GenesData)){
		GenesData$Type <- NULL
	}
}
GeneMean <- mean(GenesData[[5]])
GenesData[,"Relation"] <- GenesData[,"Mean"]/GeneMean

for (i in 1:nrow(GenesData)){
	if (GenesData[i,"Relation"] > 3){
		GenesData[i,"Type"] <- "Dup"
	} else if (GenesData[i,"Relation"] <= 0.1){
		GenesData[i,"Type"] <- "Del"
	} else {
		GenesData[i,"Type"] <- "Normal"
	}
}

Name <- str_replace(args[1],"_GenesRD.txt","_Genes_types.txt") 
write_tsv(GenesData,Name,col_names=TRUE)

