
setwd("~/Desktop/IBD/BGI_analysis")
require("DESeq2")
library("dplyr")
library("pheatmap")
library("RColorBrewer")
library("vsn")
library("ggplot2")

source("R/deseq_functions.R")
theme_set(theme_minimal() )  # pre-set the bw theme.

data<-get_data("data/gene_count_matrix_mouse_annot.csv", "data/config_mouse.txt")
cts<-data[[1]]; coldata<-data[[2]]

data<-get_data_miRNA("data/miRNAs_expressed_all_samples_1536160042.csv", "data/config_mouse.txt", 24)
#data<-get_data_miRNA("data/miRNAs_expressed_all_samples_02_07_2018_t_15_45_27.csv", "data/config_mouse.txt", 24)
cts<-data[[1]]; coldata<-data[[2]]

#blood_samples=c("M1-3","M2-3","M4-3","M7-3","M8-3","M9-3")
dds <- get_deseq(cts, coldata, c(19:24), 1)
vsd <- vst(dds, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(dds, vsd, "mouse", "colon")
resOrdered_uc_co <- get_results(dds, "mouse", "colon", 0.1)
# PCA plot
org="mouse"; tissue="colon"; type="small"
plotPCA(vsd, intgroup=c("condition")) + 
  aes(label = name,size=9)+
  geom_text_repel()+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),plot.title = element_text(hjust = 0.5,size=20,face="bold"))
ggsave(paste("figures/PCA_", org, "_",tissue,"_",type, ".png", sep=""), width=9, height=7, dpi = 100)
plotCounts(dds, "MSTRG.19359")


#plotCounts(dds, "MSTRG.19359")
blood_samples=rownames(coldata[(coldata$condition=="C" | coldata$time=="3") & coldata$tissue=="blood",])
dds <- get_deseq(cts, coldata, blood_samples, 1)
#vsd <- vst(dds, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(dds, vsd, "mouse", "blood")
resOrdered_uc_co <- get_results(dds, "mouse", "blood", 0.1)
org="mouse"; tissue="blood";type="small"
plotPCA(vsd, intgroup=c("condition")) + 
  aes(label = name,size=9)+
  geom_text_repel()+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),plot.title = element_text(hjust = 0.5,size=20,face="bold"))
ggsave(paste("figures/PCA_", org, "_",tissue,"_",type, ".png", sep=""), width=9, height=7, dpi = 100)
#plotCounts(dds,"ENSMUSG00000072875")


source("R/deseq_functions.R")

data<-get_data("data/gene_count_matrix_pig_annot.csv", "data/config_pig.txt")
cts<-data[[1]]; coldata<-data[[2]]

data<-get_data_miRNA("data/miRNAs_expressed_all_samples_1536167449.csv", "data/config_pig.txt", 30)
cts<-data[[1]]; coldata<-data[[2]]

#colon_samples=rownames(coldata[coldata$tissue=="colon" & rownames(coldata)!="PC10" & rownames(coldata)!="PC16" ,]) 
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
dds <- get_deseq(cts, coldata, colon_samples, 1)
#vsd <- vst(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(dds, vsd, "pig", "colon")
resOrdered_uc_co <- get_results(dds, "pig", "colon", 0.05)
org="pig"; tissue="colon";type="small"
plotPCA(vsd, intgroup=c("condition")) + 
  aes(label = name,size=9)+
  geom_text_repel()+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),plot.title = element_text(hjust = 0.5,size=20,face="bold"))
ggsave(paste("figures/PCA_", org, "_",tissue,"_",type, ".png", sep=""), width=9, height=7, dpi = 100)


blood_samples=rownames(coldata[(coldata$condition=="C" | coldata$time=="4") & coldata$tissue=="blood",])
#blood_samples=rownames(coldata[coldata$time=="4" & coldata$tissue=="blood",])
#blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="D",])

dds <- get_deseq(cts, coldata, blood_samples, 1)
vsd <- vst(dds, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
meanSdPlot(assay(vsd))
basic_plots(dds, vsd, "pig", "blood")
resOrdered_uc_co <- get_results(dds, "pig", "blood", 0.05)
org="pig"; tissue="blood";
plotPCA(vsd, intgroup=c("condition")) + 
#plotPCA(vsd, intgroup=c("time")) + 
  aes(label = name,size=9)+
  geom_text_repel()+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),plot.title = element_text(hjust = 0.5,size=20,face="bold"))
ggsave(paste("figures/PCA_", org, "_",tissue,".png", sep=""), width=9, height=7, dpi = 100)






