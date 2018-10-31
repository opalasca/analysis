setwd("~/Desktop/IBD/BGI_analysis")

source("R/deseq_functions.R")
#theme_set(theme_minimal() )  # pre-set the bw theme.

m_orth_file="data/mouse_pig_human_orth.tsv"
p_orth_file="data/pig_human_mouse_orth.tsv"
m_orth <- read.csv(m_orth_file, header=TRUE, sep="\t")
p_orth <- read.csv(p_orth_file, header=TRUE, sep="\t")
mh_orth <- get_orth_table(m_orth, 5, c(1,3,4)) 
mp_orth <- get_orth_table(m_orth, 8, c(1,6,7))
ph_orth <- get_orth_table(p_orth, 5, c(1,3,4)) 
pm_orth <- get_orth_table(p_orth, 9, c(1,7,8))

data<-get_data("data/gene_count_matrix_mouse_annot.csv", "data/config_mouse.txt")
ctsm<-data[[1]]; coldatam<-data[[2]]
coldatam$batch <- rep("mouse",24)
ctsm <- cbind(rownames(ctsm), data.frame(ctsm, row.names=NULL,check.names=FALSE))
colnames(ctsm)[1]<-"id"

data<-get_data("data/gene_count_matrix_pig_annot.csv", "data/config_pig.txt")
#data<-get_data("data/pig_quant_fc.csv", "data/config_pig.txt")
ctsp<-data[[1]]; coldatap<-data[[2]]
coldatap$batch <- rep("pig",30)
ctsp <- cbind(rownames(ctsp), data.frame(ctsp, row.names=NULL,check.names=FALSE))
colnames(ctsp)[1]<-"id"

ctsm_orth <- as.data.frame(merge(ctsm, mp_orth[,c(1,2)], by.x='id', by.y='Gene.stable.ID'))
cts_merged <- as.data.frame(merge(ctsm_orth, ctsp, by.x='Pig.gene.stable.ID', by.y='id'))
rownames(cts_merged) <- cts_merged[,1]
cts_merged <- cts_merged[ -c(1:2)]
coldata_merged <- rbind(coldatam, coldatap)

colon_samples=rownames(coldata_merged[coldata_merged$tissue=="colon",]) 
ddsc <- get_deseq_mp(cts_merged, coldata_merged, colon_samples, 3)
vsd <- vst(ddsc, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
meanSdPlot(assay(vsd))

blood_samples=rownames(coldata_merged[coldata_merged$tissue=="blood" & coldata_merged$condition=="D" & 
                                        coldata_merged$batch=="pig" &
                                        (coldata_merged$time=="1" | coldata_merged$time=="3" | coldata_merged$time=="4"),])
ddsb <- get_deseq_blood_paired(cts_merged, coldata_merged, blood_samples, 1)
vsd <- vst(ddsb, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
meanSdPlot(assay(vsd))


library("PerformanceAnalytics")
my_data <- ccs[, c(1,2,3,7,8,9)]
my_data <- as.matrix(log(counts(ddsc)+1))
pdf(paste("figures/", "correlation_mouse_pig_colon", ".pdf", sep=""))
chart.Correlation(my_data, histogram=TRUE, pch=19)
dev.off()

my_data <- as.matrix(log(counts(ddsb)+1))
pdf(paste("figures/", "correlation_pig_blood_t1_t3_t4", ".pdf", sep=""))
chart.Correlation(my_data, histogram=TRUE, pch=19)
dev.off()


library(Hmisc)
ccs <- as.matrix(assay(vsd))
ccs <- as.matrix()
p<-rcorr(ccs, type="pearson")
s<-rcorr(ccs,type="spearman")
