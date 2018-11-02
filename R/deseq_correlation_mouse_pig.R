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


#Transform avg intensities per group to percentiles 
avg_intensities <- as.data.frame(fit.unique[["coefficients"]])
avg_intensities$mean <- fit.unique$Amean 
avg_intensities <- merge(avg_intensities,annot["Gene.Symbol"],by="row.names",all.x=TRUE)

data<-get_data("data/pig_strtie_quant_tpm.csv", "data/config_pig.txt")
tpmp<-data[[1]]; coldatap<-data[[2]]
data<-get_data("data/mouse_strtie_quant_tpm.csv", "data/config_mouse.txt")
tpmm<-data[[1]]; coldatam<-data[[2]]

tpmm <- transform_tpm(tpmm, coldatam, "mouse") 
tpmp <- transform_tpm(tpmp, coldatap, "pig") 

getDE <- function(tpm, samp1, samp2){
  #tpm<-tpmp;samp1="P20-1"; samp2=
  res <- tpm[,c("id", samp1, samp2, "Human.gene.stable.ID", "Human.gene.name")]
  res$logFC <- as.numeric(tpm$samp2) / as.numeric(tpm$samp1)
  return(res)
} 

resP20 <- getDE(tpmp,"P20-1","P20-4")
#avg_intensities <- as.data.frame(fit.unique[["coefficients"]])
#avg_intensities$mean <- fit.unique$Amean 

IPA <- read.csv('data/IBD_IPA.txt', sep="\t")
IBD_Reza <- read.csv('data/IBD_genes_Reza.txt', sep="\t", header=FALSE)
IBD_GWAS <- read.csv('data/IBD_associated_loci_Reza.txt', sep="\t", header=FALSE)
common_hm <- data.frame()
common_hp <- data.frame()
common_mp <- data.frame()
#res1<-as.data.frame(merge(avg_intensities, AI_vs_C, by.x='Gene.Symbol', by.y='Gene.Symbol'))
res1<-avg_intensities
res1<-as.data.frame(merge(avg_intensities, IPA, by.x='Gene.Symbol', by.y='Molecule.Name'))
res2<-tpmm
res3<-tpmp
common_hm <- as.data.frame(merge(res1, res2, by.x='Gene.Symbol', by.y='Human.gene.name'))
common_hp <- as.data.frame(merge(res1, res3, by.x='Gene.Symbol', by.y='Human.gene.name'))
common_mp <- as.data.frame(merge(res2, res3, by.x='Human.gene.name', by.y='Human.gene.name'))
common_mp_filtered <- as.data.frame(merge(common_mp, resm_filtered, by.x='Human.gene.name', by.y='Human.gene.name'))

get_correlations(common_hm, "human", "mouse") 
get_correlations(common_hp, "human", "pig") 

common_mp = common_mp_filtered
cor(common_mp$qCD.x, common_mp$qCD.y)
cor(common_mp$qCD.x, common_mp$qCC.y)
cor(common_mp$qCC.x, common_mp$qCC.y)
cor(common_mp$qCC.x, common_mp$qCD.y)

require("lsa")
cosine(as.matrix(common_hm[,c("AI","C","CD","CC")]))
cor(as.matrix(common_hm[,c("AI","C","CD","CC")]))
cosine(as.matrix(common_hp[,c("AI","C","CD","CC")]))
cor(as.matrix(common_hp[,c("AI","C","CD","CC")]))
cosine(as.matrix(common_mp[,c("CD.x","CC.x","CD.y","CC.y")]))
cor(as.matrix(common_mp[,c("CD.x","CC.x","CD.y","CC.y")]))

cosine(as.matrix(common_mp[,c("BD.x","BC.x","BD.y","BC.y")]))
cor(as.matrix(common_mp[,c("BD.x","BC.x","BD.y","BC.y")]))

drops <- c("Human.gene.stable.ID.x","Human.gene.stable.ID.y","Human.gene.name","id.x","id.y")
cos_mouse_pig <- cosine(as.matrix(common_mp[, !(names(common_mp) %in% drops)]))
write.table(cos_mouse_pig, file="results/mouse_pig_cosine.tsv", row.names=T,col.names=NA, sep="\t")




png(paste("figures/hist_AI.png", sep=""))
qplot(AI, data=avg_intensities, geom='histogram')
dev.off()
png(paste("figures/hist_C.png", sep=""))
qplot(C, data=avg_intensities, geom='histogram')
dev.off()

png(paste("figures/hist_avg_col_mouse.png", sep=""))
qplot(CD, data=tpmm, geom='histogram')
dev.off()
png(paste("figures/hist_avg_col_pig.png", sep=""))
qplot(CD, data=tpmp, geom='histogram')
dev.off()


