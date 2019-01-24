setwd("~/Desktop/IBD/BGI_analysis")
source("R/deseq_functions.R")
require(gdata)
library(data.table)
require("lsa")

mmu_name_dict = read.csv("data/coexpression/mouse_ens_refseq_entrez.tsv", header=TRUE, sep="\t")
mmu_ens_refseq <- mmu_name_dict[mmu_name_dict$RefSeq.mRNA.ID != "",][,c(1,3)]
mmu_ens_refseq <- unique(mmu_ens_refseq)
mmu_ens_entrez <- mmu_name_dict[!is.na(mmu_name_dict$NCBI.gene.ID),][,c(1,4)]
mmu_ens_entrez <- unique(mmu_ens_entrez)

mmu_ens_gene_prot <- read.csv("data/ensembl_gene_protein.tsv", header=TRUE, sep="\t")
mmu_ens_gene_prot <- mmu_ens_gene_prot[mmu_ens_gene_prot$Protein.stable.ID!="",]
mmu_ens_gene_prot <- unique(mmu_ens_gene_prot[,c(1,3)])

m_seq <- read.csv("data/mir_name_sequence_mus_musculus.csv", header=FALSE, sep="\t")
p_seq <- read.csv("data/mir_name_sequence_sus_scrofa.csv", header=FALSE, sep="\t")
colnames(m_seq)=c("id","seq")
colnames(p_seq)=c("id","seq")

# Rain interactions
mmu_rain_combined <- fread("data/coexpression/10090.v1.combined.tsv", header=FALSE, sep="\t")
mmu_rain_combined <- mmu_rain_combined[,c(2:5)]
colnames(mmu_rain_combined) <- c("miRNA", "prot", "type", "score")
mmu_rain_combined$miRNA <- tolower(mmu_rain_combined$miRNA)

mmu_rain <- merge(mmu_rain_combined, mmu_ens_gene_prot, by.x="prot", by.y="Protein.stable.ID")
mmu_rain <- mmu_rain[order(mmu_rain$Gene.stable.ID, mmu_rain$miRNA, -mmu_rain$score),]
mmu_rain <- mmu_rain[!duplicated(mmu_rain[,c("miRNA","Gene.stable.ID")]),]

#Read expression/correlation data
org="mouse"
resg<-resmb; resg2<-resmb2; resmir<-resmbs;resmir2<-resmbs2; rlog_small<-rlog_mouse_blood_small; rlog_total<-rlog_mouse_blood_total; exprg<-meantpmb ; exprmir<-meanrldbs;tissue="blood"
resg<-resmc; resg2<-NULL; resmir<-resmcs; resmir2<-NULL; rlog_small<-rlog_mouse_colon_small; rlog_total<-rlog_mouse_colon_total;  exprg<-meantpmc; exprmir<-meanrldcs; tissue="colon"

rlog_small<-rlog_small[rownames(rlog_small) %in% resmir[resmir$avg.expr>30,]$id,]
rlog_total<-rlog_total[rownames(rlog_total) %in% resg[resg$avg.expr>30,]$id,]

correl_all <- get_gene_correlations(rlog_small,rlog_total)

de_mirs <- return_sig_no_biotype(resmir, NULL, NULL, 0, 0.05, 1)
de_mirs_up <-de_mirs[de_mirs$logFC > 0.1,]
de_mirs_down <-de_mirs[de_mirs$logFC < - 0.1,]

correl_DE <- correl_all[correl_all$miRNA %in% de_mirs$id, ]
corr_targets_DE <- merge(correl_DE, mmu_rain, by.x=c("gene","miRNA"),by.y=c("Gene.stable.ID","miRNA"), all.x=TRUE)

pdf(paste("figures/hist_correl_vs_interactions_RAIN_",org,"_",tissue,".pdf",sep=""))
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE), 
       widths=c(1,1,1), heights=c(1,1,1))
hist(corr_targets_DE[is.na(corr_targets_DE$score),]$corr, main = paste("DE miRNAs and ","all genes, no targets",sep=""),
     xlab = "Spearman's correlation coefficients", xlim = range(-1:1))
hist(corr_targets_DE[corr_targets_DE$score>=0.0001,]$corr, main = paste("DE miRNAs and ","all genes, targets",sep=""),
     xlab = "Spearman's correlation coefficients", xlim = range(-1:1))
hist(corr_targets_DE[corr_targets_DE$score>=0.1,]$corr, main = paste("DE miRNAs and ","all genes, score > 0.1",sep=""),
     xlab = "Spearman's correlation coefficients", xlim = range(-1:1))
hist(corr_targets_DE[corr_targets_DE$score>=0.2,]$corr, main = paste("DE miRNAs and ","all genes, score > 0.2",sep=""),
     xlab = "Spearman's correlation coefficients", xlim = range(-1:1))
hist(corr_targets_DE[corr_targets_DE$score>=0.5,]$corr, main = paste("DE miRNAs and ","all genes, score > 0.5",sep=""),
     xlab = "Spearman's correlation coefficients", xlim = range(-1:1))
hist(corr_targets_DE[corr_targets_DE$score>=0.7,]$corr, main = paste("DE miRNAs and ","all genes, score > 0.7",sep=""),
     xlab = "Spearman's correlation coefficients", xlim = range(-1:1))
dev.off()

target_counts <- corr_targets_DE[!is.na(corr_targets_DE$score),][,c("miRNA")]
target_counts <- aggregate(list(no_targets=rep(1,nrow(target_counts))), target_counts, length)

