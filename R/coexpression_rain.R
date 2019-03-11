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

#RIsearch2 near perfect complementarity interactions
mmu_risearch <- fread("data/coexpression/risearch_mmu.tsv", header=FALSE, sep="\t")
names(mmu_risearch) = c("mirna","startpos","gene","transcript","energy","cigar")

#Read expression/correlation data
org="mouse"
resg<-resmb; resg2<-resmb2; resmir<-resmbs;resmir2<-resmbs2; rlog_small<-rlog_mouse_blood_small; rlog_total<-rlog_mouse_blood_total; exprmir<-meanrldbs;tissue="blood";#exprg<-meantpmb ; 
resg<-resmc; resg2<-NULL; resmir<-resmcs; resmir2<-NULL; rlog_small<-rlog_mouse_colon_small; rlog_total<-rlog_mouse_colon_total; exprmir<-meanrldcs; tissue="colon"; #exprg<-meantpmc;

rlog_small<-rlog_small[rownames(rlog_small) %in% resmir[resmir$avg.expr>30,]$id,]
rlog_total<-rlog_total[rownames(rlog_total) %in% resg[resg$avg.expr>30,]$id,]

de_mirs <- return_sig_no_biotype(resmir, NULL, NULL, 0, 0.2, 1)
correl_all <- get_gene_correlations(rlog_small,rlog_total)
correl_risearch <- merge(correl_all, mmu_risearch, by.x=c("gene","miRNA"),by.y=c("gene","mirna"))
correl_risearch_DE <- correl_risearch[correl_risearch$miRNA %in% de_mirs$id, ]
hist(correl_risearch_DE$corr)

logFCthr=1
de_genes <- return_sig_no_biotype(resg, NULL, NULL, 0, 0.1, 0.1)
de_genes_up <-de_genes[de_genes$logFC > logFCthr,]
de_genes_down <-de_genes[de_genes$logFC < -logFCthr,]
de_mirs <- return_sig_no_biotype(resmir, NULL, NULL, 0, 0.1, 0.4)
de_mirs_up <-de_mirs[de_mirs$logFC > logFCthr,]
de_mirs_down <-de_mirs[de_mirs$logFC < - logFCthr,]

pred_targets <- unique(mmu_rain[,c("Gene.stable.ID")]) 
pred_mirs <- unique(mmu_rain[,c("miRNA")])

# filter correlation table for genes common with pred file
correl_filtered <- correl_all[correl_all$gene %in% pred_targets$Gene.stable.ID  & correl_all$miRNA %in% pred_mirs$miRNA, ]
correl_targets <- merge(correl_filtered, mmu_rain, by.x=c("gene","miRNA"),by.y=c("Gene.stable.ID","miRNA"), all.x=TRUE)

corr_targets_rain_only <- merge(correl_filtered, mmu_rain, by.x=c("gene","miRNA"),by.y=c("Gene.stable.ID","miRNA"))
resgf<-resg[,c("id","logFC","padj")]
resmirf<-resmir[,c("id","logFC","padj")]
corr_targets_rain_genelogFC <- merge(corr_targets_rain_only, resgf, by.x="gene", by.y="id")
colnames(corr_targets_rain_genelogFC)[c(8,9)] <- c("gene_logFC","gene_padj")
corr_targets_rain_logFC <- merge(corr_targets_rain_genelogFC, resmirf, by.x="miRNA", by.y="id")
colnames(corr_targets_rain_logFC)[c(10,11)] <- c("mir_logFC","mir_padj")
write.table(corr_targets_rain_logFC, file=paste("results/corr_targets_rain_logFC_", org, "_", tissue, ".tsv", sep=""), sep="\t", row.names = FALSE, quote = F)


correl_DE <- correl_filtered[correl_filtered$miRNA %in% de_mirs$id & correl_filtered$gene %in% de_genes$id, ]
correl_targets_DE <- merge(correl_DE, mmu_rain, by.x=c("gene","miRNA"),by.y=c("Gene.stable.ID","miRNA"), all.x=TRUE)

pred <- dim(correl_targets_DE[!is.na(correl_targets_DE$score),])[1]

pred_same_way <- dim(correl_targets_DE[!is.na(correl_targets_DE$score) & 
        ((correl_targets_DE$miRNA %in% de_mirs_up$id & correl_targets_DE$gene %in% de_genes_up$id) | 
         (correl_targets_DE$miRNA %in% de_mirs_down$id & correl_targets_DE$gene %in% de_genes_down$id)),] )[1]
pred_opposite <- dim(correl_targets_DE[!is.na(correl_targets_DE$score) & 
        ((correl_targets_DE$miRNA %in% de_mirs_up$id & correl_targets_DE$gene %in% de_genes_down$id) | 
          (correl_targets_DE$miRNA %in% de_mirs_down$id & correl_targets_DE$gene %in% de_genes_up$id)),] )[1]
not_pred_same_way <- dim(correl_targets_DE[is.na(correl_targets_DE$score) & 
                                         ((correl_targets_DE$miRNA %in% de_mirs_up$id & correl_targets_DE$gene %in% de_genes_up$id) | 
                                            (correl_targets_DE$miRNA %in% de_mirs_down$id & correl_targets_DE$gene %in% de_genes_down$id)),] )[1]
not_pred_opposite <- dim(correl_targets_DE[is.na(correl_targets_DE$score) & 
                                         ((correl_targets_DE$miRNA %in% de_mirs_up$id & correl_targets_DE$gene %in% de_genes_down$id) | 
                                            (correl_targets_DE$miRNA %in% de_mirs_down$id & correl_targets_DE$gene %in% de_genes_up$id)),] )[1]

m <- matrix(c(pred_same_way,pred_opposite,not_pred_same_way,not_pred_opposite), nrow=2)
colnames(m) <- c("pred target","not pred target")
rownames(m) <- c("same way","opposite")
fisher.test(m)
kable(m)





m <- matrix(c(1000,1100,10000,12000), nrow=2)
(pred_same_way/pred_opposite)/(not_pred_same_way/not_pred_opposite)



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

target_counts <- corr_targets_DE[corr_targets_DE$score>0.1,][,c("score","miRNA")]
target_counts <- target_counts[!is.na(target_counts$score),][,c("miRNA")]
target_counts <- aggregate(list(no_targets=rep(1,nrow(target_counts))), target_counts, length)

#TO DO:
# Produce same plot as above with the pairs weighted by expression differences between DSS and Ctrl for 
# miRNAs and genes 

target_counts_expr <- merge(target_counts, exprmir, by.x="miRNA",by.y="id")
corr_targets_mir_down_avg_expr <- merge(corr_targets_mir_down, exprmir, by.x="miRNA",by.y="id")

hist(target_counts$no_targets)  


# Take the set1 of up or down regulated mirs and compute the predicted ranked targets of set1  
correl_mir_up <- correl_all[correl_all$miRNA %in% de_mirs_up$id, ]
correl_mir_down <- correl_all[correl_all$miRNA %in% de_mirs_down$id, ]

# DE miRNAs and all totalRNA genes found in our data, with target info from RAIN (not all pairs have target relationship)
corr_targets_mir_up <- merge(correl_mir_up, mmu_rain, by.x=c("gene","miRNA"),by.y=c("Gene.stable.ID","miRNA"))
corr_targets_mir_down <- merge(correl_mir_down, mmu_rain, by.x=c("gene","miRNA"),by.y=c("Gene.stable.ID","miRNA"))

corr_targets_mir_up_avg_expr <- merge(corr_targets_mir_up, exprmir, by.x="miRNA",by.y="id")
corr_targets_mir_down_avg_expr <- merge(corr_targets_mir_down, exprmir, by.x="miRNA",by.y="id")

targets_mir_up_avg_expr <- corr_targets_mir_up_avg_expr[,list(score=sum(score*Diff)), by=c("gene")]
targets_mir_down_avg_expr <- corr_targets_mir_down_avg_expr[,list(score=sum(-score*Diff)), by=c("gene")]

expr_targets_mir_up <- merge(targets_mir_up_avg_expr, de_genes, by.x='gene',by.y='id',all.x=TRUE,all.y=TRUE)
expr_targets_mir_down <- merge(targets_mir_down_avg_expr, de_genes, by.x='gene',by.y='id',all.x=TRUE,all.y=TRUE)

expr_targets_mir_up[is.na(expr_targets_mir_up$score),]$score<-0
expr_targets_mir_down[is.na(expr_targets_mir_down$score),]$score<-0
expr_targets_combined <- (merge(expr_targets_mir_up, expr_targets_mir_down, by='gene',all.x=TRUE,all.y=TRUE))
expr_targets_combined$diff_up_down <- expr_targets_combined$score.x - expr_targets_combined$score.y
expr_targets_combined <- expr_targets_combined[expr_targets_combined$diff_up_down!=0,]

pairs <- expr_targets_combined[!is.na()]

pdf(paste("figures/target_loc_vs_logFC_combined","_",org,"_",tissue,".pdf",sep=""))
#layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), widths=c(1,1), heights=c(1,1))
plot(jitter(expr_targets_combined$diff_up_down), expr_targets_combined$logFC.x, pch=16,cex=0.4,main="Targets of DE miRNAs",ylab="target logFC",xlab="Targeting score, s.up-s.down")
#plot(jitter(expr_targets_combined$diff_up_down), expr_targets_combined$logFC.y, pch=16,cex=0.4,main="Targets of up-reg miRNAs",ylab="logFC",xlab="Target score")
dev.off()

pdf(paste("figures/target_loc_vs_logFC","_",org,"_",tissue,".pdf",sep=""))
layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), widths=c(1,1), heights=c(1,1))
plot(jitter(expr_targets_mir_up$score), expr_targets_mir_up$logFC, pch=16,cex=0.4,main="Targets of up-reg miRNAs",ylab="logFC",xlab="Target score")
plot(jitter(expr_targets_mir_down$score), expr_targets_mir_down$logFC, pch=16,cex=0.4,main="Targets of down-reg miRNAs",ylab="logFC",xlab="Target score")
dev.off()
