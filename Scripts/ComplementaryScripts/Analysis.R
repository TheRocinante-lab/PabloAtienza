
library("readr")
library(stringr)
#Extraction of data
args=commandArgs(trailingOnly=TRUE)

Data <- read_tsv(args[1], skip = 0, col_names = TRUE)
if ("Size" %in% colnames(Data)){
Size <- Data[,"Size", drop = FALSE] 
} else {
Size <- Data[,"End",drop =FALSE]- Data[,"Start",drop = FALSE] +1
}
Ovlp <- Data[,"Number of bases overlapping RLCRs", drop = FALSE]
Percentage <- (Ovlp/Size)*100 ; Data[, "Percentage of overlapping bases"] <- Percentage

SvtyData <- Data[c(Data[,"Percentage of overlapping bases"] < 85),]

SvtyName <- str_replace(args[1],".cnv.RLCR","_85_filtered.txt") 

write_tsv(Data,args[1], col_names= TRUE)
write_tsv(SvtyData,SvtyName,col_names= TRUE)
