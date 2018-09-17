
setwd("~/Desktop/IBD/BGI_analysis")
require("DESeq2")
library("dplyr")
library("pheatmap")
library("RColorBrewer")
library("vsn")
library("ggplot2")

source("R/deseq_functions.R")
theme_set(theme_minimal() )  # pre-set the bw theme.

m_orth_file="data/mouse_pig_human_orth.tsv"
p_orth_file="data/pig_human_mouse_orth.tsv"
m_orth <- read.csv(m_orth_file, header=TRUE, sep="\t")
p_orth <- read.csv(p_orth_file, header=TRUE, sep="\t")

mh_orth <- get_orth_table(m_orth, 5, c(1,3,4)) 
mp_orth <- get_orth_table(m_orth, 8, c(1,6,7))
ph_orth <- get_orth_table(p_orth, 5, c(1,3,4)) 
pm_orth <- get_orth_table(p_orth, 9, c(1,7,8))

data<-get_data("data/gene_count_matrix_mouse_annot.csv", "data/config_mouse.txt")
cts<-data[[1]]; coldata<-data[[2]]
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
dds <- get_deseq(cts, coldata, colon_samples, 3)
vsd <- vst(dds, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
head(assay(vsd), 3)
meanSdPlot(assay(vsd))
basic_plots(dds, vsd, "mouse", "colon","total")
resOrdered_uc_co <- get_results(dds, "mouse", "colon", "total", 0.1)
plot_gene("ENSMUSG00000072875",dds)


data<-get_data("data/gene_count_matrix_mouse_annot.csv", "data/config_mouse.txt")
cts<-data[[1]]; coldata<-data[[2]]
blood_samples=rownames(coldata[(coldata$condition=="C" | coldata$time=="3") & coldata$tissue=="blood",])
#blood_samples=rownames(coldata[coldata$tissue=="blood",])
dds <- get_deseq(cts, coldata, blood_samples, 3)
vsd <- vst(dds, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
meanSdPlot(assay(vsd))
basic_plots(dds, vsd, "mouse", "blood", "total")
resOrdered_uc_co <- get_results(dds, "mouse", "blood","total", 0.1)


plotCounts(dds,"ENSMUSG00000072875")
plot_gene("ENSMUSG00000072875",dds)


topGene="ENSMUSG00000072875"
d <- plotCounts(dds, gene=topGene, intgroup="condition", returnData=TRUE)
# Plotting the gene normalized counts, using the samplenames (rownames of d as labels)
ggplot(d, aes(x = condition, y = count, color = condition)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle(topGene) +
  theme(plot.title = element_text(hjust = 0.5))

#p <- plotCounts(dds, gene=topGene, aes(label=colnames(dds)))
ggplot(pdata, aes(label=colnames(dds))) 

ggplot(data, aes(x=condition, y=count, fill=condition)) 
 # scale_y_log10()  

#scale_y_log10() 
# geom_point(position=position_jitter(width=.1,height=0))


data<-get_data_miRNA("data/miRNAs_expressed_all_samples_1536160042.csv", "data/config_mouse.txt", 24)
#data<-get_data_miRNA("data/miRNAs_expressed_all_samples_02_07_2018_t_15_45_27.csv", "data/config_mouse.txt", 24)
cts<-data[[1]]; coldata<-data[[2]]
dds <- get_deseq(cts, coldata, c(19:24), 3)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(dds, vsd, "mouse", "colon","small")
resOrdered_uc_co <- get_results(dds, "mouse", "colon", "small", 0.1)

data<-get_data_miRNA("data/miRNAs_expressed_all_samples_1536160042.csv", "data/config_mouse.txt", 24)
#data<-get_data_miRNA("data/miRNAs_expressed_all_samples_02_07_2018_t_15_45_27.csv", "data/config_mouse.txt", 24)
cts<-data[[1]]; coldata<-data[[2]]
blood_samples=rownames(coldata[(coldata$condition=="C" | coldata$time=="3") & coldata$tissue=="blood",])
dds <- get_deseq(cts, coldata, blood_samples, 3)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(dds, vsd, "mouse", "blood", "small")
resOrdered_uc_co <- get_results(dds, "mouse", "blood","small", 0.1)

source("R/deseq_functions.R")

data<-get_data("data/gene_count_matrix_pig_annot.csv", "data/config_pig.txt")
cts<-data[[1]]; coldata<-data[[2]]
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
dds <- get_deseq(cts, coldata, colon_samples, 3)
#vsd <- vst(dds, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
vsd <- vst(dds, blind=TRUE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
meanSdPlot(assay(vsd))
basic_plots(dds, vsd, "pig", "colon","total")
resOrdered_uc_co <- get_results(dds, "pig", "colon", "total", 0.1)

data<-get_data("data/gene_count_matrix_pig_annot.csv", "data/config_pig.txt")
cts<-data[[1]]; coldata<-data[[2]]
#blood_samples=rownames(coldata[(coldata$condition=="C" | coldata$time=="4") & coldata$tissue=="blood",])
blood_samples=rownames(coldata[coldata$time=="4" & coldata$tissue=="blood",])
blood_samples=rownames(coldata[coldata$tissue=="blood",])
#blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="D",])
dds <- get_deseq(cts, coldata, blood_samples, 3)
vsd <- vst(dds, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
meanSdPlot(assay(vsd))
basic_plots(dds, vsd, "pig", "blood", "total")
resOrdered_uc_co <- get_results(dds, "pig", "blood","total", 0.1)


data<-get_data_miRNA("data/miRNAs_expressed_all_samples_1536167449.csv", "data/config_pig.txt", 30)
cts<-data[[1]]; coldata<-data[[2]]
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
dds <- get_deseq(cts, coldata, colon_samples, 3)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(dds, vsd, "pig", "colon","small")
resOrdered_uc_co <- get_results(dds, "pig", "colon", "small", 0.1)

blood_samples=rownames(coldata[(coldata$condition=="C" | coldata$time=="4") & coldata$tissue=="blood",])
dds <- get_deseq(cts, coldata, blood_samples, 3)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(dds, vsd, "pig", "blood","small")
resOrdered_uc_co <- get_results(dds, "pig", "blood", "small", 0.1)

#plotCounts(dds, "MSTRG.19359")



