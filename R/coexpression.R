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

m_seq <- read.csv("data/mir_name_sequence_mus_musculus.csv", header=FALSE, sep="\t")
p_seq <- read.csv("data/mir_name_sequence_sus_scrofa.csv", header=FALSE, sep="\t")
colnames(m_seq)=c("id","seq")
colnames(p_seq)=c("id","seq")

# MiRNA host genes, computed by intersecting mirbase gff files with ensembl gff files
mmu_host = read.csv("data/coexpression/mir_host_mouse.tsv", header=FALSE, sep="\t")
#ssc_host = read.csv("data/coexpression/ssc_host_mouse.tsv", header=FALSE, sep="\t")

# Validated targets from mirTarBase
# Entrez id 
#mmu_validated = read.xls("data/coexpression/mmu_MTI.xls", sheet = 1, header = TRUE)
#ssc_validated = read.xls("data/coexpression/ssc_MTI.xls", sheet = 1, header = TRUE)
#mmu_validated$miRNA <- tolower(mmu_validated$miRNA)

# Predicted targets from miRWalk
mmu_predicted <- fread("data/coexpression/mmu_miRWalk_3UTR.txt", header=TRUE, sep="\t")
#mmu_pred <- mmu_predicted[,1:3]
mmu_pred <- merge(mmu_predicted[,1:3], mmu_ens_refseq, by.x="mRNA", by.y="RefSeq.mRNA.ID")
mmu_pred$miRNA <- tolower(mmu_pred$miRNA)

cos.sim=function(ma, mb){
  mat=tcrossprod(ma, mb)
  t1=sqrt(apply(ma, 1, crossprod))
  t2=sqrt(apply(mb, 1, crossprod))
  mat / outer(t1,t2)
}

get_gene_correlations <- function(df1,df2,mmu_pred){
  #df1 small rna, fd2 long rna
  df1<-df1[!grepl("^.+(mt|st)$",rownames(df1)),]
  c<-(cor(t(df2), t(df1), method="spearman"))
  #c<-(cos.sim(df2, df1))
  c<-melt(c)
  c<-data.table(c, keep.rownames = TRUE)
  names(c)<-c("rn","gene","miRNA","corr")
  #fwrite(c, "data/coexpression/corr_mouse_colon.csv", sep="\t")
  return(c)  
}

# Mirna targets coexpression
org="mouse"
pred_full<-mmu_pred
pred<-aggregate(list(targetloc=rep(1,nrow(pred_full))), pred_full, length)

resg<-resmb; resg2<-resmb2; resmir<-resmbs;resmir2<-resmbs2; rlog_small<-rlog_mouse_blood_small; rlog_total<-rlog_mouse_blood_total; tissue="blood"
resg<-resmc; resg2<-NULL; resmir<-resmcs; resmir2<-NULL; rlog_small<-rlog_mouse_colon_small; rlog_total<-rlog_mouse_colon_total; tissue="colon"

correl_all <- get_gene_correlations(rlog_small,rlog_total)
pred_targets <- unique(pred[,c("Gene.stable.ID")]) 
pred_mirs <- unique(pred[,c("miRNA")])
# restrict the sets of genes and miRNAs to those found in miRWalk predictions
common_mirs <- merge(unique(pred_full[,c("miRNA")]), unique(correl_all[,c("miRNA")]), by="miRNA")
common_genes <- merge(unique(pred_full[,c("Gene.stable.ID")]), unique(correl_all[,c("gene")]), by.x="Gene.stable.ID",by.y="gene")
# filter correlation table for genes common with pred file
correl <- correl_all[correl_all$gene %in% common_genes$Gene.stable.ID & correl_all$miRNA %in% common_mirs$miRNA, ]
#filter correlation table for DE genes and mirs
corr_targets <- merge(correl,pred,by.x=c("gene","miRNA"),by.y=c("Gene.stable.ID","miRNA"), all.x=TRUE)

de_genes <- return_sig_no_biotype(resg, resg2, NULL, 0, 0.1, 0.1)
de_genes_up <-de_genes[de_genes$logFC > 1,]
de_genes_down <-de_genes[de_genes$logFC < -1,]
de_mirs <- return_sig_no_biotype(resmir, resmir2, NULL, 0, 0.05, 1)
de_mirs_up <-de_mirs[de_mirs$logFC > 0.1,]
de_mirs_down <-de_mirs[de_mirs$logFC < - 0.1,]

deg<-de_genes_down; dem<-de_mirs_up; dirg="down"; dirmir="up"; 
#deg<-de_genes_up; dem<-de_mirs_down; dirg="up"; dirmir="down"; 

correl_DE_sp <- correl[correl$gene %in% deg$id & correl$miRNA %in% dem$id, ]
correl_DE <- correl[correl$gene %in% de_genes$id & correl$miRNA %in% de_mirs$id, ]
corr_targets_DE <- merge(correl_DE,pred,by.x=c("gene","miRNA"),by.y=c("Gene.stable.ID","miRNA"), all.x=TRUE)
corr_targets_DE[is.na(corr_targets_DE$targetloc),]$targetloc<-0
plot(jitter(corr_targets_DE$corr), jitter(corr_targets_DE$targetloc), pch=16,cex=0.4)

pdf(paste("figures/hist_miRNA_targets_deg_targetloc",dirg,"_mir_",dirmir,"_",org,"_",tissue,".pdf",sep=""))
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE), 
       widths=c(1,1,1), heights=c(1,1,1))
hist(corr_targets_DE[corr_targets_DE$targetloc==0,]$corr, main = paste("DE miRNAs and ","DE genes, no targets",sep=""),
     xlab = "Spearman's correlation coefficients", xlim = range(-1:1))
hist(corr_targets_DE[corr_targets_DE$targetloc==1,]$corr, main = paste("DE miRNAs and ","DE genes, 1 target",sep=""),
     xlab = "Spearman's correlation coefficients", xlim = range(-1:1))
hist(corr_targets_DE[corr_targets_DE$targetloc==2,]$corr, main = paste("DE miRNAs and ","DE genes, 2 targets",sep=""),
     xlab = "Spearman's correlation coefficients", xlim = range(-1:1))
hist(corr_targets_DE[corr_targets_DE$targetloc==3,]$corr, main = paste("DE miRNAs and ","DE genes, 3 targets",sep=""),
     xlab = "Spearman's correlation coefficients", xlim = range(-1:1))
hist(corr_targets_DE[corr_targets_DE$targetloc==4,]$corr, main = paste("DE miRNAs and ","DE genes, 4 targets",sep=""),
     xlab = "Spearman's correlation coefficients", xlim = range(-1:1))
hist(corr_targets_DE[corr_targets_DE$targetloc==5,]$corr, main = paste("DE miRNAs and ","DE genes, 5 targets",sep=""),
     xlab = "Spearman's correlation coefficients", xlim = range(-1:1))
dev.off()

# Take the set1 of up or down regulated mirs and compute the predicted ranked targets of set1  
correl_mir_up <- correl[correl$miRNA %in% de_mirs_up$id, ]
correl_mir_down <- correl[correl$miRNA %in% de_mirs_down$id, ]

# DE miRNAs and all totalRNA genes found in our data, with target info (not all pairs have target relationship)
corr_targets_mir_up_full <- merge(correl_mir_up, pred, by.x=c("gene","miRNA"),by.y=c("Gene.stable.ID","miRNA"), all.x=TRUE)
corr_targets_mir_down_full <- merge(correl_mir_down, pred, by.x=c("gene","miRNA"),by.y=c("Gene.stable.ID","miRNA"), all.x=TRUE)

# DE miRNAs and their targets, also found in our data 
corr_targets_mir_up <- merge(correl_mir_up, pred, by.x=c("gene","miRNA"),by.y=c("Gene.stable.ID","miRNA"))
corr_targets_mir_down <- merge(correl_mir_down, pred, by.x=c("gene","miRNA"),by.y=c("Gene.stable.ID","miRNA"))

# get a set of ranked putative targets
targets_mir_up <- corr_targets_mir_up[,list(targetloc=sum(targetloc)), by=c("gene","Genesymbol")]
targets_mir_down <- corr_targets_mir_down[,list(targetloc=sum(targetloc)), by=c("gene","Genesymbol")]
expr_targets_mir_up <- merge(targets_mir_up, de_genes, by.x='gene',by.y='id',all.x=TRUE,all.y=TRUE)
expr_targets_mir_down <- merge(targets_mir_down, de_genes, by.x='gene',by.y='id',all.x=TRUE,all.y=TRUE)

plot(jitter(expr_targets_mir_up$targetloc), expr_targets_mir_up$logFC, pch=16,cex=0.4)
plot(jitter(expr_targets_mir_down$targetloc), expr_targets_mir_down$logFC, pch=16,cex=0.4)

# how many of the genes from the putative ranked target set are down/up-reg

dim(corr_targets_mir_up)[1]
length(de_genes_down$id)
dim(corr_targets_mir_up[corr_targets_mir_up$gene %in% de_genes_down$id,])[1]
dim(corr_targets_mir_up[corr_targets_mir_up$gene %in% de_genes_up$id,])[1]

# Out of the up-reg mirs, how many of the down-reg genes are targets?
dim(corr_targets_mir_up[corr_targets_mir_up$gene %in% de_genes_down$id & 
                          !is.na(corr_targets_mir_up$Genesymbol),])[1]/
  dim(corr_targets_mir_up[corr_targets_mir_up$gene %in% de_genes_down$id,])[1]

# Out of the up-reg mirs, how many of the up-reg genes are targets?
dim(corr_targets_mir_up[corr_targets_mir_up$gene %in% de_genes_up$id & 
                          !is.na(corr_targets_mir_up$Genesymbol),])[1]/
  dim(corr_targets_mir_up[corr_targets_mir_up$gene %in% de_genes_up$id,])[1]

# Out of the down-reg mirs, how many of the down-reg genes are targets?
dim(corr_targets_mir_down[corr_targets_mir_down$gene %in% de_genes_down$id & 
                            !is.na(corr_targets_mir_down$Genesymbol),])[1]/
  dim(corr_targets_mir_down[corr_targets_mir_down$gene %in% de_genes_down$id,])[1]

# how many of the up-reg genes are targets of the down-reg mirs, ?
dim(corr_targets_mir_down[corr_targets_mir_down$gene %in% de_genes_up$id & 
                            !is.na(corr_targets_mir_down$Genesymbol),])[1]/
  dim(corr_targets_mir_down[corr_targets_mir_down$gene %in% de_genes_up$id,])[1]


corrt <- corr_targets_DE
#corrt <- corr_targets
pdf(paste("figures/hist_miRNA_targets_deg",dirg,"_mir_",dirmir,"_",org,"_",tissue,".pdf",sep=""))
layout(matrix(c(1,1,2,2,3,3), 3, 2, byrow = TRUE), 
       widths=c(1,1,1), heights=c(1,1,1))
hist(corrt$corr, main = paste(dirmir,"-regulated miRNAs and ",dirg, "-regulated genes",sep=""),
     xlab = "Pearson's correlation coefficients", xlim = range(-1:1))
hist(corrt$corr[!is.na(corrt$Genesymbol)],
     main = paste(dirmir,"-regulated miRNAs and ",dirg, "-regulated TARGET genes",sep=""),
     xlab = "Pearson's correlation coefficients", xlim = range(-1:1))
hist(corrt$corr[is.na(corrt$Genesymbol)],
     main = paste(dirmir,"-regulated miRNAs and ",dirg, "-regulated non-TARGET genes",sep=""),
     xlab = "Pearson's correlation coefficients", xlim = range(-1:1))
dev.off()

corr_targets <- merge(correl,pred,by.x=c("gene","miRNA"),by.y=c("Gene.stable.ID","miRNA"))

length(corrt$corr[!is.na(corrt$Genesymbol)])/length(corrt$corr)
length(corrt$corr[!is.na(corrt$Genesymbol)])/length(corr_targets$corr)

target_vs_nontarget_ratio
# miRNA-host coexpression
# Precursors are used in host table -> add a column to the host matrix with -5p and -3p included
# (if they exist in some mature miRNA list?)

host <- mmu_host
names(host)<-c("miRNA","gene","gname","biotype") 
host_5p <- host
host_3p <- host
host_5p$mature <- paste(host$miRNA, "-5p", sep="")
host_3p$mature <- paste(host$miRNA, "-3p", sep="")
host <- rbind(host_5p,host_3p)
host <- merge(host, m_seq, by.x="mature", by.y="id")

corr_host <- merge(correl_all, host, by.x=c("gene","miRNA"),by.y=c("gene","mature"))
corr_host_targets <- merge(corr_targets, host, by.x=c("gene","miRNA"),by.y=c("gene","mature"))
corr_host_DE <- merge(correl_DE, host, by.x=c("gene","miRNA"),by.y=c("gene","mature"))

pdf(paste("figures/hist_miRNA_host_targets_", org,"_",tissue,".pdf",sep=""))
layout(matrix(c(1,1,2,2,3,3), 3, 2, byrow = TRUE), 
       widths=c(1,1,1), heights=c(1,1,1))
hist(corr_host$corr, main = "miRNAs and their host genes", xlab = "Spearman's correlation coefficients",xlim = range(-1:1))
hist(corr_host_DE$corr, main = "DE miRNAs and their DE host genes", xlab = "Spearman's correlation coefficients",xlim = range(-1:1))
hist(corr_host_targets$corr, main =  "miRNAs and their host genes which are also predicted as targets", xlab = "Spearman's correlation coefficients",xlim = range(-1:1))
#hist(corr_host_targets$corr[!is.na(corr_host_targets$Genesymbol)], main = "miRNAs and their host genes, which occur in the predicted target list", xlab = "Pearson's correlation coefficients",xlim = range(-1:1))
dev.off()














# Commands to obtain mirna-host pairs from mirbase download files
mirna_host_mouse = read.csv("data/mirna_context_mouse.txt", header=FALSE, sep="\t")
mirna_host_pig = read.csv("data/mirna_context_pig.txt", header=FALSE, sep="\t")
prec_id_name <- read.csv("data/mirna.txt", header=FALSE, sep="\t")
gene_tr_name_mouse <- read.csv("data/mouse_gene_transcript_name.tsv", header=TRUE, sep="\t")
gene_tr_name_pig <- read.csv("data/mouse_gene_transcript_name.tsv", header=TRUE, sep="\t")
mirna_host <- mirna_host_mouse[,c(1:2,4)]
names(mirna_host)<- c("internal_id","host_tr_id","region")
id_name <- prec_id_name[,c(1,3,4,5)]
names(id_name)<- c("internal_id","name","alternative_name","full_name")
m_host <- as.data.frame(merge(mirna_host, id_name, by="internal_id",all.x=TRUE ))

corr_targets_mir_up_f <- merge(correl_mir_up, pred, by.x=c("gene","miRNA"),by.y=c("Gene.stable.ID","miRNA"))
corr_targets_mir_down_f <- merge(correl_mir_down, pred, by.x=c("gene","miRNA"),by.y=c("Gene.stable.ID","miRNA"))

dim(corr_targets_mir_up_f[corr_targets_mir_up_f$gene %in% de_genes_down$id,])[1]/ dim(corr_targets_mir_up_f)[1]
dim(corr_targets_mir_up_f[corr_targets_mir_up_f$gene %in% de_genes_up$id,])[1]/ dim(corr_targets_mir_up_f)[1]

dim(corr_targets_mir_down_f[corr_targets_mir_down_f$gene %in% de_genes_up$id,])[1]/ 
  dim(corr_targets_mir_down[corr_targets_mir_down$gene %in% de_genes_up$id,])[1]
dim(corr_targets_mir_down_f[corr_targets_mir_down_f$gene %in% de_genes_down$id,])[1]/ 
  dim(corr_targets_mir_down[corr_targets_mir_down$gene %in% de_genes_down$id,])[1]

length(unique(corr_targets_mir_down_f[corr_targets_mir_down_f$gene %in% de_genes_up$id,]$gene))/length(de_genes_up$id)
length(unique(corr_targets_mir_up_f[corr_targets_mir_up_f$gene %in% de_genes_up$id,]$gene))/length(de_genes_up$id)

length(unique(corr_targets_mir_down_f[corr_targets_mir_down_f$gene %in% de_genes_down$id,]$gene))/length(de_genes_down$id)
length(unique(corr_targets_mir_up_f[corr_targets_mir_up_f$gene %in% de_genes_down$id,]$gene))/length(de_genes_down$id)


length((corr_targets_mir_down_f[corr_targets_mir_down_f$gene %in% de_genes_up$id,]$gene))/dim(corr_targets_mir_down[corr_targets_mir_down$gene %in% de_genes_up$id,])[1]

length(unique(corr_targets_mir_down_f[corr_targets_mir_down_f$gene %in% de_genes_down$id,]$gene))/length(corr_targets_mir_down_f$gene)

length(unique(corr_targets_mir_up_f[corr_targets_mir_up_f$gene %in% de_genes_up$id,]$gene))/length(corr_targets_mir_up_f$gene)
length(unique(corr_targets_mir_up_f[corr_targets_mir_up_f$gene %in% de_genes_down$id,]$gene))/length(corr_targets_mir_up_f$gene)

