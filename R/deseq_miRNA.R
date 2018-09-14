#install.packages("Vennerable", repos="http://R-Forge.R-project.org") 
#source("https://bioconductor.org/biocLite.R")
#biocLite("affy")
#biocLite("DESeq2")
#biocLite("pheatmap")

setwd("~/Desktop/IBD/BGI_analysis")
require("DESeq2")
library("dplyr")

cts <- read.csv(file=paste("data/",acc,"_miRNAs_expressed_all_samples_26_04_2018_t_13_43_32.csv",sep=""), sep="\t")
#remove duplicate mature miRNAs 
cts<-cts[order(cts$X.miRNA,-cts$read_count),]
cts<-cts[!duplicated(cts$X.miRNA),]
#cts <- as.matrix(cts)
cts<-cts[,c(1,5:41)]
# Assign the first column as row names
cts2 <- cts[,-1]
cts2 <- as.matrix(cts2)
#rownames(cts2) <- cts[,1]
rownames(cts2) <- tolower(gsub("_.*","",cts[,1]))

coldata <- read.csv(file=paste("data/",acc,"_config2.txt",sep=""), sep="\t", header=FALSE)
coldata <- coldata[,c(2,3,4)]
coldata2 <- coldata[,-1]
rownames(coldata2) <- coldata[,1]
colnames(coldata2) <- c("condition","type")

#arrange columns of cts in the same order as rows of colnames 
all(rownames(coldata2) %in% colnames(cts2))
cts2 <- cts2[, rownames(coldata2)]
all(rownames(coldata2) == colnames(cts2))

dds <- DESeqDataSetFromMatrix(countData = cts2,
                              colData = coldata2,
                              design = ~ condition)
dds <- DESeq(dds)

#thr=1

res_UC_DD <- results(dds, contrast=c("condition","UC","DD"))
res_CD_DD <- results(dds, contrast=c("condition","CD","DD"))
res_UC_CD <- results(dds, contrast=c("condition","UC","CD"))

resOrdered_UC_DD <- as.data.frame(res_UC_DD[order(res_UC_DD$pvalue),])
resOrdered_CD_DD <- as.data.frame(res_CD_DD[order(res_CD_DD$pvalue),])
resOrdered_UC_CD <- as.data.frame(res_UC_CD[order(res_UC_CD$pvalue),])

names(resOrdered_UC_DD)[names(resOrdered_UC_DD) == "log2FoldChange"] = "logFC"
names(resOrdered_CD_DD)[names(resOrdered_CD_DD) == "log2FoldChange"] = "logFC"
names(resOrdered_UC_CD)[names(resOrdered_UC_CD) == "log2FoldChange"] = "logFC"

names(resOrdered_UC_DD)[names(resOrdered_UC_DD) == "baseMean"] = "avg.expr"
names(resOrdered_CD_DD)[names(resOrdered_CD_DD) == "baseMean"] = "avg.expr"
names(resOrdered_UC_CD)[names(resOrdered_UC_CD) == "baseMean"] = "avg.expr"

#resOrdered_UC_DD$mirna <-rownames(resOrdered_UC_DD)
#resOrdered_CD_DD$mirna <-rownames(resOrdered_CD_DD)
#resOrdered_UC_CD$mirna <-rownames(resOrdered_UC_CD)

resOrdered_UC_DD <- subset(resOrdered_UC_DD, padj < thr)
resOrdered_CD_DD <- subset(resOrdered_CD_DD, padj < thr)
resOrdered_UC_CD <- subset(resOrdered_UC_CD, padj < thr)

resOrdered_UC_DD$common_id <-tolower(gsub("_.*","",rownames(resOrdered_UC_DD)))
resOrdered_CD_DD$common_id <-tolower(gsub("_.*","",rownames(resOrdered_CD_DD)))
resOrdered_UC_CD$common_id <-tolower(gsub("_.*","",rownames(resOrdered_UC_CD)))

write.table(resOrdered_UC_DD[c(1,2,6,7)], file=paste("results/",acc,"_UC_vs_DD_p_", thr, ".tsv", sep=""), sep="\t")
write.table(resOrdered_CD_DD[c(1,2,6,7)], file=paste("results/",acc,"_CD_vs_DD_p_", thr, ".tsv", sep=""), sep="\t")
write.table(resOrdered_UC_CD[c(1,2,6,7)], file=paste("results/",acc,"_UC_vs_CD_p_", thr, ".tsv",sep=""), sep="\t")

#common_1_2 <- as.data.frame(merge(resOrdered_UC_DD[c(1,2,6,7)], resOrdered_CD_DD[c(1,2,6,7)], by=c('mirna','mirna'))[c(1,2,3,4,6,7)])
#common_1_3 <- as.data.frame(merge(resOrdered_UC_DD[c(1,2,6,7)], resOrdered_UC_CD[c(1,2,6,7)], by=c('mirna','mirna'))[c(1,2,3,4,6,7)])
#common_2_3 <- as.data.frame(merge(resOrdered_CD_DD[c(1,2,6,7)], resOrdered_UC_CD[c(1,2,6,7)],  by=c('mirna','mirna'))[c(1,2,3,4,6,7)])

