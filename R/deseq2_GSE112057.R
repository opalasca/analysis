setwd("~/Desktop/IBD/BGI_analysis")
source("R/deseq_functions.R")

counts_file = "data/human_GSE112057/GSE112057_RawCounts_dataset.txt"
config_file = "data/human_GSE112057/config_GSE112057.txt"
cts <- read.csv(counts_file, check.names=FALSE,header=TRUE, sep="\t")
#cts<-cts[order(cts$Gene),]
rownames(cts) <- make.names(cts$Gene, unique = TRUE)
cts2 <- cts
cts2 <- cts2[,-1]
rownames(cts2) <- make.names(cts$Gene, unique = TRUE)
cts <- as.matrix(cts2)

conditions_to_keep=c("CD","UC","Control")
coldata <- read.csv(config_file, sep="\t")
coldata$condition <- as.character(coldata$condition)
coldata$condition[coldata$condition == "Crohn's Disease"] <- "CD"
coldata$condition[coldata$condition == "Ulcerative Colitis"] <- "UC"
coldata<-as.data.frame(coldata)
coldata <- coldata[coldata$condition %in% conditions_to_keep,]
rownames(coldata)=coldata$sample
coldata2=as.data.frame(coldata[,c("condition")])
rownames(coldata2)=coldata$sample
colnames(coldata2)=c("condition")
coldata2$condition <- as.factor(coldata2$condition)
cts<-cts[,colnames(cts) %in% rownames(coldata)]
cts <- cts[, rownames(coldata)]
print(all(rownames(coldata) %in% colnames(cts)))

sample_list=colnames(cts)
sample_list=rownames(coldata)
rowsums=1
thr=0
ddsh <- DESeqDataSetFromMatrix(countData = cts[,sample_list],
                                colData = coldata[sample_list,],
                                design = ~ condition)
print(nrow(ddsh))
ddsh <- estimateSizeFactors(ddsh)
mnc <- rowMeans(counts(ddsh, normalized=TRUE))
ddsh <- ddsh[mnc > 3,]
print(nrow(ddsh))
ddsh$condition <- relevel(ddsh$condition, ref = "Control")
ddsh <- DESeq(ddsh)

vsd <- vst(ddsh, blind=TRUE) ##
rld <- rlog(ddsh, blind=TRUE)

plotPCA(vsd, intgroup=c("condition")) 
  
basic_plots(ddsh, vsd, "human", "blood","total",FALSE)
resUC <- results(ddsh, independentFiltering = TRUE, test="Wald",contrast=c("condition","UC","Control"))
resCD <- results(ddsh, independentFiltering = TRUE, test="Wald",contrast=c("condition","CD","Control"))
resUCvsCD <- results(ddsh, independentFiltering = TRUE, test="Wald",contrast=c("condition","UC","CD"))

resUC <- as.data.frame(resUC)
resCD <- as.data.frame(resCD)
resUCvsCD <- as.data.frame(resUC)

resUC$name <- rownames(resUC)
resCD$name <- rownames(resCD)
resUCvsCD$name <- rownames(resUCvsCD)
names(resUC)[names(resUC) == "log2FoldChange"] = "logFC"
names(resUC)[names(resUC) == "baseMean"] = "avg.expr"
names(resCD)[names(resCD) == "log2FoldChange"] = "logFC"
names(resCD)[names(resCD) == "baseMean"] = "avg.expr"
names(resUCvsCD)[names(resUCvsCD) == "log2FoldChange"] = "logFC"
names(resUCvsCD)[names(resUCvsCD) == "baseMean"] = "avg.expr"

sigUC = resUC[abs(resUC$logFC)>1.5 & resUC$padj<0.05,]
sigCD = resCD[abs(resCD$logFC)>1 & resCD$padj<0.05,]

sigUCvsCD = resUCvsCD[abs(resUCvsCD$logFC)>1.5 & resUCvsCD$padj<0.05,]

draw_volcano("human","blood", "total", resUC, 0.05, 1.5,"")

ddsmc<-ddsm
rlog_mouse_colon_total <- assay(rldc)
tpm_mouse_colon<-tpm_mouse[rownames(tpm_mouse) %in% rownames(rlog_mouse_colon_total),]
meanrldc<-get_mean(rlog_mouse_colon_total,coldata,"mouse", "colon", "total")
#meantpmc<-get_mean(tpm_mouse_colon, coldata, "mouse", "colon", "total")
#clust<-heatmap_DE(sig,rldc,"mouse","colon","total_protein_coding","Protein coding genes",2)
#go_clust1 <- go_term_enrichment_list(clust,resmc,0.1,"mouse","colon","1")
#go_clust2 <- go_term_enrichment_list(clust,resmc,0.1,"mouse","colon","2")
sig<-return_sig(resmc, NULL, NULL, 1, 0.05, 0.05, protein_coding)
rld<-rldc