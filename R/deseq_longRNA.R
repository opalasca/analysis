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

biotypes_to_remove <- c("misc_RNA","miscRNA", "Mt_rRNA","rRNA")
mouse_genes_to_keep <- m_biotypes[!(m_biotypes$biotype %in% biotypes_to_remove),]
pig_genes_to_keep <- p_biotypes[!(p_biotypes$biotype %in% biotypes_to_remove),]

lncRNAs=c("lincRNA","antisense", "3prime_overlapping_ncRNA","bidirectional_promoter_lncRNA",
          "processed_transcript", "sense_intronic","sense_overlapping")

data<-get_data("data/gene_count_matrix_mouse_annot.csv", "data/config_mouse.txt")
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
dim(cts)
cts<-cts[rownames(cts) %in% mouse_genes_to_keep$id,]
dim(cts)
colon_samples=rownames(coldata[coldata$tissue=="colon",])
#ddsm <- get_deseq_ref_D(cts, coldata, colon_samples, 3)
ddsm<-NULL
ddsm <- get_deseq(cts, coldata, colon_samples, 1, 0)
vsd <- vst(ddsm, blind=FALSE) ##
rldc <- rlog(ddsm, blind=TRUE)
basic_plots(ddsm, rldc, "mouse", "colon","total",FALSE)
resmc<-NULL
resmc <- as.data.frame(get_results(ddsm, "mouse", "colon", "total", 1, mh_orth, mp_orth, m_biotypes, NULL,""))
ddsmc<-ddsm
res<-resmc

lfcthr=1; pthr=0.1
siglncc <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype %in% lncRNAs,]
sigpcc <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype=="protein_coding",]
heatmap_DE(sigpcc,rldc,"mouse","colon","total_protein_coding","Protein coding genes")
heatmap_DE(siglncc,rldc,"mouse","colon","total_lncRNA","LncRNAs")
dim(sigpcc[sigpcc$logFC>1,])
dim(sigpcc[sigpcc$logFC<1,])
dim(siglncc[siglncc$logFC<1,])
dim(siglncc[siglncc$logFC>1,])

#Mouse blood day 8 DSS vs controls
data<-get_data("data/gene_count_matrix_mouse_annot_reseq.csv", "data/config_mouse.txt")
#data<-get_data("data/gene_count_matrix_mouse_annot.csv", "data/config_mouse.txt")
cts<-data[[1]]; coldata<-data[[2]]
cts<-cts[rownames(cts) %in% mouse_genes_to_keep$id,]
blood_samples=rownames(coldata[(coldata$time=="8") & coldata$tissue=="blood",])
ddsm <- get_deseq(cts, coldata, blood_samples, 1, 0)
vsd <- vst(ddsm, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
rldb <- rlog(ddsm, blind=TRUE)
meanSdPlot(assay(vsd))
basic_plots(ddsm, rldb, "mouse", "blood", "total",FALSE)
resmb<-NULL
resmb <- get_results(ddsm, "mouse", "blood","total", 1, mh_orth, mp_orth, m_biotypes,NULL,"")
ddsmb <- ddsm
boxplot(assay(rldb))
res<-resmb
lfcthr=1; pthr=0.1
siglncb <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype %in% lncRNAs,]
sigpcb <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype=="protein_coding",]
heatmap_DE(sigpcb,rldb,"mouse","blood","total_protein_coding","Protein coding genes")
heatmap_DE(siglncb,rldb,"mouse","blood","total_lncRNA","LncRNAs")
dim(sigpcb[sigpcb$logFC>1,])
dim(sigpcb[sigpcb$logFC<1,])
dim(siglncb[siglncb$logFC<1,])
dim(siglncb[siglncb$logFC>1,])

#Mouse blood time series, DSS treated
data<-get_data("data/gene_count_matrix_mouse_annot_reseq.csv", "data/config_mouse.txt")
cts<-data[[1]]; coldata<-data[[2]]
cts<-cts[rownames(cts) %in% mouse_genes_to_keep$id,]
blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" ,])
ddsm <- get_deseq_blood_time(cts, coldata, blood_samples, 1)
vsd <- vst(ddsm, blind=FALSE) ##
rldb <- rlog(ddsm, blind=TRUE); #ntd <- normTransform(dds)
meanSdPlot(assay(vsd))
basic_plots(ddsm, rldb, "mouse", "blood", "total",TRUE)
resmb<-NULL
resmb8 <- get_results(ddsm, "mouse", "blood","total_ts8vs0", 1, mh_orth, mp_orth, m_biotypes,NULL,"8")
resmb2<-NULL
resmb2 <- get_results(ddsm, "mouse", "blood","total_ts2vs0", 1, mh_orth, mp_orth, m_biotypes,NULL,"2")
ddsmb <- ddsm
boxplot(assay(rldb))
lfcthr=1; pthr=0.1
res<-resmb8
siglnc <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype %in% lncRNAs,]
sigpc <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype=="protein_coding",]
heatmap_DE_ts(sigpc,rldb,"mouse","blood","total_protein_coding","Protein coding")
heatmap_DE_ts(siglnc,rldb,"mouse","blood","total_lncRNA","LncRNAs")


data<-get_data("data/gene_count_matrix_pig_annot.csv", "data/config_pig.txt")
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
dim(cts)
cts<-cts[rownames(cts) %in% pig_genes_to_keep$id,]
dim(cts)
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
#colon_samples=rownames(coldata[coldata$tissue=="colon" & coldata$subject!="P20" & coldata$subject!="P17",]) 
#colon_samples=rownames(coldata[coldata$tissue=="colon" & coldata$subject!="P17",]) 
#colon_samples=rownames(coldata[coldata$tissue=="colon" & coldata$subject!="P20",]) 
#ddsp <- get_deseq_ref_D(cts, coldata, colon_samples, 3)
ddsp<-NULL
ddsp <- get_deseq_batch(cts, coldata, colon_samples, 1, 0)
#ddsp <- get_deseq(cts, coldata, colon_samples, 2, 0)
vsd <- vst(ddsp, blind=FALSE) 
rldc <- rlog(ddsp, blind=TRUE); 
meanSdPlot(assay(vsd))
basic_plots(ddsp, rldc, "pig", "colon","total",FALSE)
respc<-NULL
respc <- as.data.frame(get_results(ddsp, "pig", "colon", "total", 1, ph_orth, pm_orth, p_biotypes,NULL,""))
ddspc<-ddsp
boxplot(assay(rldc))
res<-respc
lfcthr=0; pthr=0.1
siglncc <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype %in% lncRNAs,]
sigpcc <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype=="protein_coding",]
heatmap_DE(sigpcc,rldc,"pig","colon","total_protein_coding","Protein coding genes")
heatmap_DE(siglncc,rldc,"pig","colon","total_lncRNA","LncRNAs")
dim(sigpcb[sigpcc$logFC>1,])
dim(sigpcb[sigpcc$logFC<1,])
dim(siglncb[siglncc$logFC>1,])
dim(siglncb[siglncc$logFC<1,])

#plotPCA3D(vsd, intgroup = "condition", ntop = 500,returnData = FALSE)

data<-get_data("data/gene_count_matrix_pig_annot.csv", "data/config_pig.txt")
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
cts<-cts[rownames(cts) %in% pig_genes_to_keep$id,]
#blood_samples=rownames(coldata[coldata$tissue=="blood",]) 
blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS",])
ddsp <- NULL
ddsp <- get_deseq_blood_time(cts, coldata, blood_samples, 1)
vsd <- vst(ddsp, blind=FALSE) 
rldb <- rlog(ddsp, blind=TRUE); 
meanSdPlot(assay(vsd))
basic_plots(ddsp, rldb, "pig", "blood", "total", TRUE)
respb<-NULL
respb4 <- as.data.frame(get_results(ddsp, "pig", "blood", "total", 1, ph_orth, pm_orth, p_biotypes,NULL,"4"))
respb2 <- as.data.frame(get_results(ddsp, "pig", "blood", "total", 1, ph_orth, pm_orth, p_biotypes,NULL,"2"))
respb5 <- as.data.frame(get_results(ddsp, "pig", "blood", "total", 1, ph_orth, pm_orth, p_biotypes,NULL,"5"))
res<-respb4
ddspb <- ddsp
boxplot(assay(rldb))
lfcthr=1; pthr=0.1
siglncb <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype %in% lncRNAs,]
sigpcb <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype=="protein_coding",]
heatmap_DE_ts(sigpcb,rldb,"pig","blood","total_protein_coding","Protein coding genes")
heatmap_DE_ts(siglncb,rldb,"pig","blood","total_lncRNA","LncRNAs")
dim(sigpcb[sigpcb$logFC>1,])
dim(sigpcb[sigpcb$logFC<1,])
dim(siglncb[siglncb$logFC>1,])
dim(siglncb[siglncb$logFC<1,])

plotPCA3D(vsd, intgroup = "subject", ntop = 500,returnData = FALSE)







