#!/usr/bin/Rscript

args <- commandArgs(TRUE)
srcFile <- args[1]
outFile <- args[2]

library(data.table)
library(sqldf)
library(stringr)

data<-fread(srcFile, header=T, na.strings=c(".", "null"))
data$TYPE <- ifelse((nchar(data$MOTHER_GT)==3 & nchar(data$FATHER_GT)==3 & nchar(data$CHILD_GT)==3), "SNP", "INDEL")

#Original MDNM settings
filtDataMdnm<-subset(data, sapply(strsplit(as.character(FATHER_PL), ","), `[`, 1) == "0" & sapply(strsplit(as.character(MOTHER_PL), ","), `[`, 1) == "0" & sapply(strsplit(as.character(CHILD_PL), ","), `[`, 2) == "0" & MOTHER_DP>=10 & FATHER_DP>=10 & CHILD_DP>=10 & MOTHER_DP<=100 & FATHER_DP<=100 & CHILD_DP<=100 & !is.na(TP) & sapply(strsplit(as.character(MOTHER_AD), ","), `[`, 2) == "0" & sapply(strsplit(as.character(FATHER_AD), ","), `[`, 2) == "0" & TP>=30 & as.numeric(sapply(strsplit(as.character(CHILD_AD), ","), `[`, 2)) >= 5 & as.numeric(sapply(strsplit(as.character(CHILD_AD), ","), `[`, 1)) >= 5 & as.numeric(sapply(strsplit(as.character(CHILD_AD), ","), `[`, 2))/CHILD_DP >= 0.3)

write.table(filtDataMdnm,paste(outFile, "mdnm.txt", sep="."),row.names=F,col.names=T,quote=F,sep="\t") 

#Epilepsy paper settings
filtDataEpi<-subset(data, MOTHER_DP>=10 & FATHER_DP>=10 & CHILD_DP>=(MOTHER_DP + FATHER_DP)/10 & as.numeric(sapply(strsplit(as.character(FATHER_AD), ","), `[`, 2))/FATHER_DP <= 0.05 & as.numeric(sapply(strsplit(as.character(MOTHER_AD), ","), `[`, 2))/MOTHER_DP <= 0.05 & as.numeric(sapply(strsplit(as.character(CHILD_AD), ","), `[`, 2))/CHILD_DP >= 0.25 & sapply(strsplit(as.character(MOTHER_PL), ","), `[`, 1) == "0" & as.numeric(sapply(strsplit(as.character(MOTHER_PL), ","), `[`, 2)) > 20 & as.numeric(sapply(strsplit(as.character(MOTHER_PL), ","), `[`, 3)) > 20 & sapply(strsplit(as.character(FATHER_PL), ","), `[`, 1) == "0" & as.numeric(sapply(strsplit(as.character(FATHER_PL), ","), `[`, 2)) > 20 & as.numeric(sapply(strsplit(as.character(FATHER_PL), ","), `[`, 3)) > 20 & sapply(strsplit(as.character(CHILD_PL), ","), `[`, 2) == "0" & as.numeric(sapply(strsplit(as.character(CHILD_PL), ","), `[`, 1)) > 20 & as.numeric(sapply(strsplit(as.character(CHILD_PL), ","), `[`, 3)) > 0 & as.numeric(sapply(strsplit(as.character(CHILD_AD), ","), `[`, 2)) >= 3) 

write.table(filtDataEpi,paste(outFile, "epi.txt", sep="."),row.names=F,col.names=T,quote=F,sep="\t")

#Relaxed MDNM settings; all except DP<=100, AD==0, TP>=30
filtDataMdnmRelaxed<-subset(data, sapply(strsplit(as.character(FATHER_PL), ","), `[`, 1) == "0" & sapply(strsplit(as.character(MOTHER_PL), ","), `[`, 1) == "0" & sapply(strsplit(as.character(CHILD_PL), ","), `[`, 2) == "0" & MOTHER_DP>=10 & FATHER_DP>=10 & CHILD_DP>=10 & as.numeric(sapply(strsplit(as.character(CHILD_AD), ","), `[`, 2)) >= 5 & as.numeric(sapply(strsplit(as.character(CHILD_AD), ","), `[`, 1)) >= 5 & as.numeric(sapply(strsplit(as.character(CHILD_AD), ","), `[`, 2))/CHILD_DP >= 0.3 & !is.na(TP))

write.table(filtDataMdnmRelaxed,paste(outFile, "mdnm.relaxed.txt", sep="."),row.names=F,col.names=T,quote=F,sep="\t")

#Get calls where C/M/F is hom-alt/ref/ref
filtHomAlt<-subset(data,sapply(strsplit(as.character(FATHER_PL), ","), `[`, 1) == "0" & sapply(strsplit(as.character(MOTHER_PL), ","), `[`, 1) == "0" & sapply(strsplit(as.character(CHILD_PL), ","), `[`, 3) == "0" & MOTHER_DP>=10 & FATHER_DP>=10 & CHILD_DP>=10)

write.table(filtHomAlt,paste(outFile, "homalt.txt", sep="."),row.names=F,col.names=T,quote=F,sep="\t")

filtDataEpi$CHROM<-as.character(filtDataEpi$CHROM)
filtDataMdnm$CHROM<-as.character(filtDataMdnm$CHROM)
common<-merge(filtDataEpi, filtDataMdnm)
write.table(common,paste(outFile, "common.txt", sep="."),row.names=F,col.names=T,quote=F,sep="\t")
