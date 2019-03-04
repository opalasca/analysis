setwd("~/Desktop/IBD/BGI_analysis")
source("R/deseq_functions.R")
library(gplots)
#theme_set(theme_minimal() )  # pre-set the bw theme.

m_orth_file="data/mouse_pig_human_orth.tsv"
p_orth_file="data/pig_human_mouse_orth.tsv"
m_orth <- read.csv(m_orth_file, header=TRUE, sep="\t")
p_orth <- read.csv(p_orth_file, header=TRUE, sep="\t")
mh_lncRNA_orth <- read.csv("data/lncRNA_mouse_human_orthologs.tsv", header=TRUE, sep="\t")
mh_lncRNA_orth <- mh_lncRNA_orth[,c(5,2,1)]
names(mh_lncRNA_orth) <- names(mh_orth)
mh_orth <- get_orth_table(m_orth, 5, c(1,3,4)) 
mh_orth <- rbind(mh_orth,mh_lncRNA_orth)
mh_orth <- unique(mh_orth)
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
protein_coding=c("protein_coding")
pc_lnc=c(lncRNAs,protein_coding)

mnames <- read.csv("data/mouse_gene_transcript_name.tsv", header=TRUE, sep="\t", stringsAsFactors = FALSE)
pnames <- read.csv("data/pig_gene_transcript_name.tsv", header=TRUE, sep="\t",stringsAsFactors = FALSE)
mnames$name <- ifelse( (is.na(m_names$Gene.name) | m_names$Gene.name==""), m_names$Gene.stable.ID, m_names$Gene.name)
mnames<-unique(mnames[c(1,6)])
pnames$name <- ifelse( (is.na(p_names$Gene.name) | p_names$Gene.name==""), p_names$Gene.stable.ID, p_names$Gene.name)
pnames<-unique(pnames[c(1,6)])
m_names <- read.csv("data/mouse_gene_transcript_name.tsv", header=TRUE, sep="\t")
p_names <- read.csv("data/pig_gene_transcript_name.tsv", header=TRUE, sep="\t")


tpmm<-get_data("data/mouse_strtie_quant_tpm_reseq.csv", "data/config_mouse.txt")
tpm_mouse<-as.data.frame(log(tpmm[[1]]+1))
hist(tpm_mouse$MC1)
tpmp<-get_data("data/pig_strtie_quant_tpm.csv", "data/config_pig.txt")
tpm_pig<-as.data.frame(log(tpmp[[1]]+1))
tpm_pig$id<-rownames(tpm_pig)
hist(tpm_pig$PC20)


data<-get_data("data/gene_count_matrix_mouse_annot.csv", "data/config_mouse.txt")
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
coldata<-coldata[order(coldata$condition,coldata$time),]
dim(cts)
cts<-cts[rownames(cts) %in% mouse_genes_to_keep$id,]
dim(cts)
colon_samples=rownames(coldata[coldata$tissue=="colon",])
#ddsm <- get_deseq_ref_D(cts, coldata, colon_samples, 3)
ddsm<-NULL
ddsm <- get_deseq(cts, coldata, colon_samples, 1, 0)
vsd <- vst(ddsm, blind=FALSE) ##
rldc <- rlog(ddsm, blind=TRUE)
#basic_plots(ddsm, rldc, "mouse", "colon","total",FALSE)
resmc<-NULL
resmc <- as.data.frame(get_results(ddsm, "mouse", "colon", "total", 1, mh_orth, mp_orth, m_biotypes, m_names, NULL,"",NULL,NULL,""))
resold <- read.csv("/Users/jmx139/Dropbox/AniGen/BGI_analysis/results/4&10&11-01-2019/mouse_colon_total_1_thr.tsv", header=TRUE, sep="\t")
resold$logFCdiff = resold$logFC - resmc$logFC
resold$pdiff = resold$padj - resmc$padj
ddsmc<-ddsm
sig<-return_sig(resmc, NULL, NULL, 1, 0.05, 0.05, protein_coding)
#heatmap_DE(sig,rldb,"mouse","blood","total_protein_coding","Protein coding genes")
#heatmap_DE(sig,rldc,"mouse","colon","total_protein_coding","Protein coding genes","")
#heatmap_DE(siglncc,rldc,"mouse","colon","total_lncRNA","LncRNAs","")
rlog_mouse_colon_total <- assay(rldc)
tpm_mouse_colon<-tpm_mouse[rownames(tpm_mouse) %in% rownames(rlog_mouse_colon_total),]
meanrldc<-get_mean(rlog_mouse_colon_total,coldata,"mouse", "colon", "total")
meantpmc<-get_mean(tpm_mouse_colon, coldata, "mouse", "colon", "total")
#clust<-heatmap_DE(sig,rldc,"mouse","colon","total_protein_coding","Protein coding genes",2)
#go_clust1 <- go_term_enrichment_list(clust,resmc,0.1,"mouse","colon","1")
#go_clust2 <- go_term_enrichment_list(clust,resmc,0.1,"mouse","colon","2")
rld<-rldc
res<-resmc
rows <- match(sig$id, row.names(rld))
mat <- assay(rld)[rows,]
mat<-t(scale(t(mat)))
heatmap_DE_kmeans(sig,rld, mat,"mouse","blood","total_protein_coding","Protein coding genes",4, FALSE)
#set.seed(42)
set.seed(7)
km<- kmeans(mat,3,iter.max = 1000) # determin how many cluster you want, I specify 2 here
m.kmeans<- cbind(mat, km$cluster) # combine the cluster with the matrix
dim(m.kmeans)[2]
o<- order(m.kmeans[,dim(m.kmeans)[2]]) # order the last column
m.kmeans<- m.kmeans[o,] # order the matrix according to the order of the last column
colnames(m.kmeans)[7] <- c("cluster")
heatmap_DE_kmeans(sig,rld,m.kmeans[,1:6],"mouse","colon","total_protein_coding","Protein coding genes",4, FALSE)
clust<-as.data.frame(m.kmeans)
clust_counts <- as.data.frame(clust[,c("cluster")])
clust_counts <- aggregate(list(cluster=rep(1,nrow(clust_counts))), clust_counts, length)
clust_counts
colnames(m.kmeans)[7] <- c("cluster")
heatmap_DE_kmeans(sig,rld,m.kmeans[,1:6],"mouse","colon","total_protein_coding","Protein coding genes",4, FALSE)
clust<-as.data.frame(m.kmeans)
clust_counts <- as.data.frame(clust[,c("cluster")])
clust_counts <- aggregate(list(cluster=rep(1,nrow(clust_counts))), clust_counts, length)
clust_counts
tclust<-clust
tclust$id <- rownames(tclust)
tclust<-tclust[,c("id","cluster")]
write.table(tclust, file=paste("results/expression_clusters_", "mouse", "_", "colon", "_", as.Date(Sys.time()), ".tsv", sep=""), sep="\t", row.names = FALSE, quote = F)
go_clust1 <- go_term_enrichment_list(clust,res,0.5,"mouse","blood","1")
go_clust2 <- go_term_enrichment_list(clust,res,0.5,"mouse","blood","2")
go_clust3 <- go_term_enrichment_list(clust,res,0.5,"mouse","blood","3")
go_clust4 <- go_term_enrichment_list(clust,res,0.5,"mouse","blood","4")
go_clust5 <- go_term_enrichment_list(clust,res,0.5,"mouse","blood","5")
sig_clust <- sig[sig$id %in% rownames(clust[clust$cluster=="1",]),]
heatmap_DE(sig_clust,rldc,"mouse","colon","total_protein_coding","Protein coding genes",2, TRUE)

# Mouse blood all samples #
###########################
data<-get_data("data/gene_count_matrix_mouse_annot_reseq.csv", "data/config_mouse.txt")
#data<-get_data("data/gene_count_matrix_mouse_annot.csv", "data/config_mouse.txt")
cts<-data[[1]]; coldata<-data[[2]]
coldata<-coldata[order(coldata$condition,coldata$time),]
cts<-cts[rownames(cts) %in% mouse_genes_to_keep$id,]
blood_samples=rownames(coldata[coldata$tissue=="blood",])
#ddsm <- get_deseq(cts, coldata, blood_samples, 1, 0)
ddsm <- get_deseq_time_cond_merged(cts, coldata, blood_samples, 1, 0,"day0_DSS")
vsd <- vst(ddsm, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
rldb <- rlog(ddsm, blind=TRUE)
#filtered_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" & coldata$time != "2",])
filtered_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" ,])
rldb <- rldb[,colnames(rldb) %in% filtered_samples]
#meanSdPlot(assay(vsd))
#basic_plots(ddsm, rldb, "mouse", "blood", "total",FALSE)
resmb<-NULL
resmb <- as.data.frame(get_results(ddsm, "mouse", "blood", "total", 1, mh_orth, mp_orth, 
                                   m_biotypes, m_names, NULL,"","day8_DSS","day0_DSS","day8DSS_vs_all_ctrl_timepoints"))
resmb2 <- as.data.frame(get_results(ddsm, "mouse", "blood", "total", 1, mh_orth, mp_orth, 
                                   m_biotypes, m_names, NULL,"","day2_DSS","day0_DSS","day2DSS_vs_all_ctrl_timepoints"))
ddsmb <- ddsm
boxplot(assay(rldb))
#sig<-return_sig(resmb,NULL, NULL, 1 , 0.05, 0.05, protein_coding)
res<-resmb[abs(resmb$logFC)>=1 & resmb$padj<=0.05,]

sig<-return_sig(resmb, resmb2, NULL, 1 , 0.05, 0.05, protein_coding)
#sig<-return_sig(resmb, NULL, NULL, 1 , 0.05, 0.1, protein_coding)
heatmap_DE_ts(sig,rldb,"mouse","blood","total_protein_coding","Protein coding genes",2)
sig<-return_sig(resmb, resmb2, NULL, 1, 0.05, 0.1, lncRNAs)
heatmap_DE_ts(sig,rldb,"mouse","blood","total_lncRNA","LncRNAs",2)
#hist(resmb$pvalue,breaks = 0:20/20)
rlog_mouse_blood_total <- assay(rldb)
tpm_mouse_blood<-tpm_mouse[rownames(tpm_mouse) %in% rownames(rlog_mouse_blood_total),]
meanrldb<-get_mean(rlog_mouse_blood_total,coldata,"mouse", "blood", "total")
meantpmb<-get_mean(tpm_mouse_blood,coldata,"mouse", "blood", "total")
#hist(meantpmb$DSS_day8)

#sig<-return_sig(resmb2, NULL, NULL, 1, 0.05, 0.5, protein_coding)
#sig<-return_sig(resmb, resmb2, NULL, 1 , 0.05, 0.1, pc_lnc)
sig<-merge(sig, mnames, by.x='id', by.y='Gene.stable.ID')
clust<-heatmap_DE_ts(sig,rldb,"mouse","blood","total_protein_coding","Protein coding genes",4)

clust_counts <- as.data.frame(clust[,c("cluster")])
clust_counts <- aggregate(list(cluster=rep(1,nrow(clust_counts))), clust_counts, length)
clust_counts
go_clust1 <- go_term_enrichment_list(clust,resmb,0.5,"mouse","blood","1")
go_clust2 <- go_term_enrichment_list(clust,resmb,0.1,"mouse","blood","2")
go_clust3 <- go_term_enrichment_list(clust,resmb,0.1,"mouse","blood","3")
go_clust4 <- go_term_enrichment_list(clust,resmb,0.1,"mouse","blood","4")
go_clust6 <- go_term_enrichment_list(clust,resmb,0.1,"mouse","blood","6")
go_clust8 <- go_term_enrichment_list(clust,resmb,0.1,"mouse","blood","8")
go_clust9 <- go_term_enrichment_list(clust,resmb,0.1,"mouse","blood","9")
go_clust11 <- go_term_enrichment_list(clust,resmb,0.1,"mouse","blood","11")
sig_clust <- sig[sig$id %in% rownames(clust[clust$cluster=="2",]),]
heatmap_DE_ts_rownames(sig_clust,rldb,"mouse","blood","total_protein_coding","Protein coding genes",2)
subclust<-heatmap_DE_ts(sig_clust,rldb,"mouse","blood","total_protein_coding","Protein coding genes",2)
go_subclust1 <- go_term_enrichment_list(subclust,resmb,0.1,"mouse","blood","1")
go_subclust2 <- go_term_enrichment_list(subclust,resmb,0.1,"mouse","blood","2")

sig_clust <- sig[sig$id %in% rownames(subclust[subclust$cluster=="2",]),]
heatmap_DE_ts_rownames(sig_clust,rldb,"mouse","blood","total_protein_coding","Protein coding genes",2)


mouse_blood_cluster_M72 <- sig[sig$id %in% rownames(clust[clust$cluster=="8",]),]

# Pig colon #
#############
data<-get_data("data/gene_count_matrix_pig_annot.csv", "data/config_pig.txt")
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
coldata<-coldata[order(coldata$condition,coldata$time),]
dim(cts)
cts<-cts[rownames(cts) %in% pig_genes_to_keep$id,]
dim(cts)
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
#ddsp <- get_deseq_ref_D(cts, coldata, colon_samples, 3)
ddsp<-NULL
ddsp <- get_deseq_batch(cts, coldata, colon_samples, 1, 0)
vsd <- vst(ddsp, blind=FALSE) 
rldc <- rlog(ddsp, blind=TRUE); 
meanSdPlot(assay(vsd))
#basic_plots(ddsp, rldc, "pig", "colon","total",FALSE)
respc<-NULL
respc <- as.data.frame(get_results(ddsp, "pig", "colon", "total", 1, ph_orth, pm_orth, p_biotypes, p_names, NULL,"",NULL,NULL,""))
ddspc<-ddsp
boxplot(assay(rldc))
sig<-return_sig(respc, NULL, NULL, 0, 0.05, 0.1, protein_coding)
sig<-merge(sig, p_names, by.x='id', by.y='Gene.stable.ID')
heatmap_DE(sig,rldc,"pig","colon","total_protein_coding","Protein coding genes",2)
#sig<-return_sig(respc, NULL, NULL, 0, 0.05, 0.1, lncRNAs)
dim(sigpcc[sigpcc$logFC>1,])
dim(sigpcc[sigpcc$logFC<1,])
dim(siglncc[siglncc$logFC>1,])
dim(siglncc[siglncc$logFC<1,])
hist(respc$pvalue)
rlog_pig_colon_total <- assay(rldc)
tpm_pig_colon<-tpm_pig[rownames(tpm_pig) %in% rownames(rlog_pig_colon_total),]
meanrldcp<-get_mean(rlog_pig_colon_total, coldata, "pig", "colon", "total")
meantpmcp<-get_mean(tpm_pig_colon, coldata, "pig", "colon", "total")
hist(meantpmcp$DSS)
clust<-heatmap_DE(sig,rldc,"pig","colon","total_protein_coding","Protein coding genes",2)
clust_counts <- as.data.frame(clust[,c("cluster")])
clust_counts <- aggregate(list(cluster=rep(1,nrow(clust_counts))), clust_counts, length)
clust_counts
#go_clust1 <- go_term_enrichment_list(clust,respb,0.1,"pig","colon","1")
#go_clust2 <- go_term_enrichment_list(clust,respb,0.1,"pig","colon","2")
go_clust3 <- go_term_enrichment_list(clust,respb,0.1,"pig","colon","1")
go_clust4 <- go_term_enrichment_list(clust,respb,0.1,"pig","colon","2")

sig_clust <- sig[sig$id %in% rownames(clust[clust$cluster=="2",]),]
heatmap_DE(sig_clust,rldc,"pig","colon","total_protein_coding","Protein coding genes",2)
pig_colon_cluster_up <- sig[sig$id %in% rownames(clust[clust$cluster=="2",]),]



#Pig blood all samples 
data<-get_data("data/gene_count_matrix_pig_annot.csv", "data/config_pig.txt")
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
coldata<-coldata[order(coldata$condition,coldata$time),]
cts<-cts[rownames(cts) %in% pig_genes_to_keep$id,]
blood_samples=rownames(coldata[coldata$tissue=="blood",]) 
ddsp <- NULL
ddsp <- get_deseq_time_cond_merged_plus_batch(cts, coldata, blood_samples, 1,0,"day0_DSS")
vsd <- vst(ddsp, blind=FALSE) 
rldb <- rlog(ddsp, blind=TRUE); 
filtered_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" ,])
rldb <- rldb[,colnames(rldb) %in% filtered_samples]
meanSdPlot(assay(vsd))
#basic_plots(ddsp, rldb, "pig", "blood", "total", TRUE)
respb<-NULL
respb <- as.data.frame(get_results(ddsp, "pig", "blood", "total", 1, ph_orth, pm_orth, 
                                   p_biotypes,p_names, NULL,"","day4_DSS","day0_DSS","day4DSS_vs_all_ctrl_timepoints"))
respb2 <- as.data.frame(get_results(ddsp, "pig", "blood", "total", 1, ph_orth, pm_orth, 
                                   p_biotypes,p_names, NULL,"","day2_DSS","day0_DSS","day2DSS_vs_all_ctrl_timepoints"))
respb5 <- as.data.frame(get_results(ddsp, "pig", "blood", "total", 1, ph_orth, pm_orth, 
                                   p_biotypes,p_names, NULL,"","day5_DSS","day0_DSS","day5DSS_vs_all_ctrl_timepoints"))
sig<-return_sig(respb, respb2, respb5, 0, 0.05, 0.1, protein_coding)
heatmap_DE_ts(sig,rldb,"pig","blood","total_protein_coding","Protein coding genes",2)
#sig<-return_sig(respb, NULL, NULL, 0, 0.05, 0.1, lncRNAs)
#heatmap_DE(sig,rldb,"pig","blood","total_lncRNA","LncRNAs","t4")
hist(respb$pvalue,breaks = 0:20/20)
rlog_pig_blood_total <- assay(rldb)
tpm_pig_blood<-tpm_pig[rownames(tpm_pig) %in% rownames(rlog_pig_blood_total),]
meanrldbp<-get_mean(rlog_pig_blood_total, coldata, "pig", "blood", "total")
meantpmbp<-get_mean(tpm_pig_blood, coldata, "pig", "blood", "total")
sig<-return_sig(respb, respb2, respb5, 1, 0.05, 0.1, pc_lnc)
sig<-merge(sig, p_names, by.x='id', by.y='Gene.stable.ID')
clust<-heatmap_DE_ts(sig,rldb,"pig","blood","total_protein_coding","Protein coding genes",8)
clust_counts <- as.data.frame(clust[,c("cluster")])
clust_counts <- aggregate(list(cluster=rep(1,nrow(clust_counts))), clust_counts, length)
clust_counts
go_clust1 <- go_term_enrichment_list(clust,respb,0.1,"pig","blood","1")
go_clust2 <- go_term_enrichment_list(clust,respb,0.1,"pig","blood","2")
go_clust3 <- go_term_enrichment_list(clust,respb,0.1,"pig","blood","3")
go_clust4 <- go_term_enrichment_list(clust,respb,0.1,"pig","blood","4")
go_clust5 <- go_term_enrichment_list(clust,respb,0.1,"pig","blood","5")
go_clust6 <- go_term_enrichment_list(clust,respb,0.1,"pig","blood","6")
go_clust7 <- go_term_enrichment_list(clust,respb,0.1,"pig","blood","7")
go_clust8 <- go_term_enrichment_list(clust,respb,0.1,"pig","blood","8")
#Draw heatmap of the desired cluster
sig_clust <- sig[sig$id %in% rownames(clust[clust$cluster=="5",]),]
heatmap_DE_ts(sig_clust,rldb,"pig","blood","total_protein_coding","Protein coding genes",2)


draw_volcano <- function(org, type, res, pthr, lfcthr){
  pthr <- 0.05 # Threshold on the adjusted p-value
  org="pig";type="total"; res=respb;lfcthr=10
  cols <- densCols(res$logFC, -log10(res$padj))
  pdf(file=paste("volcano_",org,"_",type,".pdf",sep='')) 
  plot(res$logFC, -log10(res$padj), col=cols, panel.first=grid(),
       main=paste(org,"_",type,sep=''), xlab="log2(fold-change)", ylab="-log10(adj.p-value)",
       pch=20, cex=0.6)
  abline(v=0)
  abline(v=c(-1,1), col="brown")
  abline(h=-log10(pthr), col="brown")
  selected <- abs(res$logFC) > lfcthr & res$padj < pthr 
  if (selected){
    text(res$logFC[selected],-log10(res$padj)[selected],
         lab=res$Gene.Name[selected], cex=0.4)
  }
  dev.off()
}

draw_volcano("pig","total",respb,0.01,2)

library("ggplot2") #Best plots
library("ggrepel") #Avoid overlapping labels
yourdataframe=resmb
mutateddf <- mutate(yourdataframe, sig=ifelse(yourdataframe$padj<0.1 & abs(yourdataframe$logFC)>1, "padj<0.1", "Not Sig")) #Will have different colors depending on significance
input <- cbind(gene=mutateddf$Gene.Name, mutateddf ) #convert the rownames to a column
volc = ggplot(input, aes(logFC, -log10(pvalue))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=sig)) + #add points colored by significance
  scale_color_manual(values=c("black", "red")) + 
  ggtitle("Blood") #e.g. 'Volcanoplot DESeq2'
volc+geom_text_repel(data=head(input, 50), aes(label=gene)) #adding text for the top 20 genes
ggsave("figures/Volcanoplot_blood.jpeg", device="jpeg") #In case you want to easily save to disk
volc



# Pig blood DSS samples 
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
respb4 <- as.data.frame(get_results(ddsp, "pig", "blood", "total", 1, ph_orth, pm_orth, p_biotypes,NULL,"4",NULL,NULL))
respb2 <- as.data.frame(get_results(ddsp, "pig", "blood", "total", 1, ph_orth, pm_orth, p_biotypes,NULL,"2",NULL,NULL))
respb5 <- as.data.frame(get_results(ddsp, "pig", "blood", "total", 1, ph_orth, pm_orth, p_biotypes,NULL,"5",NULL,NULL))
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

#plotPCA3D(vsd, intgroup = "subject", ntop = 500,returnData = FALSE)
#plotPCA3D(vsd, intgroup = "condition", ntop = 500,returnData = FALSE)

#Mouse blood day 8 DSS vs controls
data<-get_data("data/gene_count_matrix_mouse_annot_reseq.csv", "data/config_mouse.txt")
#data<-get_data("data/gene_count_matrix_mouse_annot.csv", "data/config_mouse.txt")
cts<-data[[1]]; coldata<-data[[2]]
coldata<-coldata[order(coldata$condition,coldata$time),]
cts<-cts[rownames(cts) %in% mouse_genes_to_keep$id,]
blood_samples=rownames(coldata[(coldata$time=="8") & coldata$tissue=="blood",])
ddsm <- get_deseq(cts, coldata, blood_samples, 1, 0)
vsd <- vst(ddsm, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
rldb <- rlog(ddsm, blind=TRUE)
meanSdPlot(assay(vsd))
resmb<-NULL
resmb <- get_results(ddsm, "mouse", "blood","total", 1, mh_orth, mp_orth, m_biotypes,m_names, NULL,"",NULL,NULL,"_day8")
ddsmb <- ddsm
#res<-resmb[abs(resmb$logFC)>=1 & resmb$padj<=0.05,]
boxplot(assay(rldb))
res<-resmb
sig<-return_sig(resmb, NULL, NULL, 1, 0.05, 0.05, protein_coding)
clust<-heatmap_DE(sig,rldb,"mouse","blood","total_protein_coding","Protein coding genes",4, TRUE)
sig<-return_sig(resmb, NULL, NULL, 1, 0.05, 0.05, lncRNAs)
heatmap_DE(sig,rldb,"mouse","blood","total_lncRNA","LncRNAs",2, TRUE)
sig<-return_sig(resmb, NULL, NULL, 1, 0.05, 0.05, protein_coding)
rld<-rldb
res<-resmb
rows <- match(sig$id, row.names(rld))
mat <- assay(rld)[rows,]
mat<-t(scale(t(mat)))
heatmap_DE_kmeans(sig,rld, mat,"mouse","blood","total_protein_coding","Protein coding genes",3, FALSE)
set.seed(42)
km<- kmeans(mat,3,iter.max=1000) # determin how many cluster you want, I specify 2 here
m.kmeans<- cbind(mat, km$cluster) # combine the cluster with the matrix
dim(m.kmeans)[2]
o<- order(m.kmeans[,dim(m.kmeans)[2]]) # order the last column
m.kmeans<- m.kmeans[o,] # order the matrix according to the order of the last column
colnames(m.kmeans)[7] <- c("cluster")
heatmap_DE_kmeans(sig,rld,m.kmeans[,1:6],"mouse","blood","total_protein_coding","Protein coding genes",4, FALSE)
clust<-as.data.frame(m.kmeans)
clust_counts <- as.data.frame(clust[,c("cluster")])
clust_counts <- aggregate(list(cluster=rep(1,nrow(clust_counts))), clust_counts, length)
clust_counts
tclust<-clust
tclust$id <- rownames(tclust)
tclust<-tclust[,c("id","cluster")]
write.table(tclust, file=paste("results/expression_clusters_", "mouse", "_", "blood", "_", as.Date(Sys.time()), ".tsv", sep=""), sep="\t", row.names = FALSE, quote = F)
go_clust1 <- go_term_enrichment_list(clust,res,0.1,"mouse","blood","1")
go_clust2 <- go_term_enrichment_list(clust,res,0.1,"mouse","blood","2")
go_clust3 <- go_term_enrichment_list(clust,res,0.1,"mouse","blood","3")
go_clust4 <- go_term_enrichment_list(clust,res,0.1,"mouse","blood","4")
sig_clust <- sig[sig$id %in% rownames(clust[clust$cluster=="3",]),]
heatmap_DE(sig_clust,rldb,"mouse","blood","total_protein_coding","Protein coding genes",2, TRUE)
#basic_plots(ddsm, rldb, "mouse", "blood", "total",FALSE)

#sig<-return_sig(resmb, NULL, NULL, 1, 0.05, 0.1, lncRNAs)
#heatmap_DE(sig,rldb,"mouse","blood","total_lncRNA","LncRNAs",2)



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
resmb8 <- get_results(ddsm, "mouse", "blood","total_ts8vs0", 1, mh_orth, mp_orth, m_biotypes,NULL,"8",NULL,NULL,"day8DSS_vs_all_ctrl_timepoints")
resmb2<-NULL
resmb2 <- get_results(ddsm, "mouse", "blood","total_ts2vs0", 1, mh_orth, mp_orth, m_biotypes,NULL,"2",NULL,NULL,"day2DSS_vs_all_ctrl_timepoints")
ddsmb <- ddsm
boxplot(assay(rldb))
lfcthr=0; pthr=0.1
res<-resmb8
siglnc <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype %in% lncRNAs,]
sigpc <- res[abs(res$logFC)>lfcthr & res$padj<pthr & res$biotype=="protein_coding",]
heatmap_DE_ts(sigpc,rldb,"mouse","blood","total_protein_coding","Protein coding")
heatmap_DE_ts(siglnc,rldb,"mouse","blood","total_lncRNA","LncRNAs")
hist(resmb8$pvalue)

