setwd("~/Desktop/IBD/BGI_analysis")
source("R/deseq_functions.R")
library(gplots)
#theme_set(theme_minimal() )  # pre-set the bw theme.

m_orth_file="data/mouse_pig_human_orth.tsv"
p_orth_file="data/pig_human_mouse_orth.tsv"
m_orth <- read.csv(m_orth_file, header=TRUE, sep="\t")
p_orth <- read.csv(p_orth_file, header=TRUE, sep="\t")
mh_orth <- get_orth_table(m_orth, 5, c(1,3,4)) 
mp_orth <- get_orth_table(m_orth, 8, c(1,6,7))
ph_orth <- get_orth_table(p_orth, 5, c(1,3,4)) 
pm_orth <- get_orth_table(p_orth, 9, c(1,7,8))

p_biotypes <- read.csv("data/biotypes/pig_biotypes.txt", header=FALSE, sep="\t")
colnames(p_biotypes) <- c("id","biotype")
m_biotypes <- read.csv("data/biotypes/mouse_biotypes.txt", header=FALSE, sep="\t")
colnames(m_biotypes) <- c("id","biotype")

data<-get_data("data/gene_count_matrix_mouse_annot.csv", "data/config_mouse.txt")
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
colon_samples=rownames(coldata[coldata$tissue=="colon",])
#ddsm <- get_deseq_ref_D(cts, coldata, colon_samples, 3)
ddsm<-NULL
ddsm <- get_deseq(cts, coldata, colon_samples, 3, 0)
vsd <- vst(ddsm, blind=FALSE) ##
rld <- rlog(ddsm, blind=FALSE)
#vsd <- rlog(ddsm, blind=FALSE);
basic_plots(ddsm, vsd, "mouse", "colon","total",FALSE)
resmc<-NULL
resmc <- as.data.frame(get_results(ddsm, "mouse", "colon", "total", 1, mh_orth, mp_orth, m_biotypes, NULL,""))
ddsmc<-ddsm
res<-resmc
lfcthr=1; pthr=0.1
siglinc <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype=="lincRNA",]
sigpc <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype=="protein_coding",]
sig <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype=="miRNA",]
heatmap_DE(sigpc,rld,"mouse","colon","total_protein_coding")
heatmap_DE(siglinc,rld,"mouse","colon","total_lincRNA")

dim(sigpc[sigpc$logFC>1,])
dim(sigpc[sigpc$logFC<1,])
dim(siglinc[siglinc$logFC<1,])
dim(siglinc[siglinc$logFC>1,])


#Mouse blood day 8 DSS vs controls
data<-get_data("data/gene_count_matrix_mouse_annot_reseq.csv", "data/config_mouse.txt")
#data<-get_data("data/gene_count_matrix_mouse_annot.csv", "data/config_mouse.txt")
cts<-data[[1]]; coldata<-data[[2]]
blood_samples=rownames(coldata[(coldata$time=="8") & coldata$tissue=="blood",])
ddsm <- get_deseq(cts, coldata, blood_samples, 3, 0)
vsd <- vst(ddsm, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
rld <- rlog(ddsm, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(ddsm, vsd, "mouse", "blood", "total",FALSE)
resmb<-NULL
resmb <- get_results(ddsm, "mouse", "blood","total", 1, mh_orth, mp_orth, m_biotypes,NULL,"")
ddsmb <- ddsm

boxplot(assay(rld))
res<-resmb
lfcthr=1; pthr=0.1
siglinc <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype=="lincRNA",]
sigpc <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype=="protein_coding",]
heatmap_DE(sigpc,rld,"mouse","blood","total_protein_coding")
heatmap_DE(siglinc,rld,"mouse","blood","total_lincRNA")

dim(sigpc[sigpc$logFC>1,])
dim(sigpc[sigpc$logFC<1,])
dim(siglinc[siglinc$logFC<1,])
dim(siglinc[siglinc$logFC>1,])

#Mouse blood time series, DSS treated
data<-get_data("data/gene_count_matrix_mouse_annot_reseq.csv", "data/config_mouse.txt")
cts<-data[[1]]; coldata<-data[[2]]
blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" ,])
ddsm <- get_deseq_blood_time(cts, coldata, blood_samples, 3)
vsd <- vst(ddsm, blind=FALSE) ##
rld <- rlog(ddsm, blind=FALSE); #ntd <- normTransform(dds)
meanSdPlot(assay(vsd))
basic_plots(ddsm, vsd, "mouse", "blood", "total",TRUE)
resmb<-NULL
resmb8 <- get_results(ddsm, "mouse", "blood","total_ts8vs0", 1, mh_orth, mp_orth, m_biotypes,NULL,"8")
resmb2<-NULL
resmb2 <- get_results(ddsm, "mouse", "blood","total_ts2vs0", 1, mh_orth, mp_orth, m_biotypes,NULL,"2")
ddsmb <- ddsm

boxplot(assay(rld))
res<-resmb8
lfcthr=1; pthr=0.1
siglinc <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype=="lincRNA",]
sigpc <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype=="protein_coding",]
heatmap_DE(sigpc,rld,"mouse","blood","total_protein_coding")
heatmap_DE(siglinc,rld,"mouse","blood","total_lincRNA")


data<-get_data("data/gene_count_matrix_pig_annot.csv", "data/config_pig.txt")
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
colon_samples=rownames(coldata[coldata$tissue=="colon" & coldata$subject!="P20" & coldata$subject!="P17",]) 
colon_samples=rownames(coldata[coldata$tissue=="colon" & coldata$subject!="P17",]) 
#colon_samples=rownames(coldata[coldata$tissue=="colon" & coldata$subject!="P20",]) 
#colon_samples=rownames(coldata[coldata$tissue=="colon" & coldata$subject!="P17",]) 
#ddsp <- get_deseq_ref_D(cts, coldata, colon_samples, 3)
ddsp<-NULL
ddsp <- get_deseq_batch(cts, coldata, colon_samples, 2, 0)
#ddsp <- get_deseq(cts, coldata, colon_samples, 2, 0)
vsd <- vst(ddsp, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
#vsd <- rlog(ddsp, blind=FALSE);
meanSdPlot(assay(vsd))
basic_plots(ddsp, vsd, "pig", "colon","total",FALSE)
respc<-NULL
respc <- as.data.frame(get_results(ddsp, "pig", "colon", "total", 1, ph_orth, pm_orth, p_biotypes,NULL,""))
ddspc<-ddsp


data<-get_data("data/gene_count_matrix_pig_annot.csv", "data/config_pig.txt")
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
#blood_samples=rownames(coldata[coldata$time=="5" & coldata$tissue=="blood",]) #use this if comparing with colon
blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" & 
                                 (coldata$time=="0" | coldata$time=="4"),])
blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS",])
ddsp <- NULL
ddsp <- get_deseq_blood_time(cts, coldata, blood_samples, 1)
vsd <- vst(ddsp, blind=FALSE) 
meanSdPlot(assay(vsd))
basic_plots(ddsp, vsd, "pig", "blood", "total", TRUE)
respb<-NULL
#respb <- results(ddsp, independentFiltering = FALSE)
respb4 <- as.data.frame(get_results(ddsp, "pig", "blood", "total", 1, ph_orth, pm_orth, p_biotypes,NULL,"4"))
respb2 <- as.data.frame(get_results(ddsp, "pig", "blood", "total", 1, ph_orth, pm_orth, p_biotypes,NULL,"2"))
respb5 <- as.data.frame(get_results(ddsp, "pig", "blood", "total", 1, ph_orth, pm_orth, p_biotypes,NULL,"5"))
ddspb <- ddsp









