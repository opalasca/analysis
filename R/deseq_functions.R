#install.packages("Vennerable", repos="http://R-Forg> install.packages("ggplot2")e.R-project.org") 
#source("https://bioconductor.org/biocLite.R")
#biocLite("hexbin")
#biocLite("vsn")
#biocLite("affy")
#biocLite("DESeq2")
#biocLite("pheatmap")
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#read-counting-step

setwd("~/Desktop/IBD/BGI_analysis")
require("DESeq2")
library("dplyr")
library("pheatmap")
library("RColorBrewer")
library("vsn")
library("ggplot2")
library(ggrepel)
library(cowplot)
#library("installr")

theme_set(theme_cowplot(font_size=12)) # reduce default font size

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Function return two dataframes - counts and config file, in a list
get_data <- function(counts_file, config_file){
  cts <- data.frame()
  coldata <- data.frame()
  cts <- as.matrix(read.csv(counts_file, row.names="gene_id", check.names=FALSE, header=TRUE))
  coldata <- read.csv(config_file, sep="\t", row.names=1)
  coldata$time <- as.factor(coldata$time)
  if ("block" %in% colnames(coldata)){
    coldata$block <- as.factor(coldata$block)
  }
  print(all(rownames(coldata) %in% colnames(cts)))
  cts <- cts[, rownames(coldata)]
  print(all(rownames(coldata) == colnames(cts)))
  data <- list(cts, coldata)
  return(data)  
}

get_data_miRNA <- function(counts_file, config_file, n){
  cts <- data.frame()
  coldata <- data.frame()
  cts <- read.csv(counts_file, check.names=FALSE, header=TRUE, sep="\t")
  colnames(cts)[1]<-"miRNA"
  cts<-cts[order(cts$miRNA,-cts$read_count),]
  cts<-cts[!duplicated(cts$miRNA),]
  cts<-cts[,c(1,5:(n+4))]
  # Assign the first column as row names
  cts2 <- cts[,-1]
  cts <- as.matrix(cts)
  rownames(cts2) <- cts[,1]
  sample_names <- read.csv("data/config_mirdeep2.txt", sep="\t", header=FALSE)
  if (n==24){
    colnames(cts2) <- sample_names[1:n,1]
  }
  else{
    colnames(cts2) <- sample_names[25:(n+24),1]
  }
  coldata <- read.csv(config_file, sep="\t", row.names=1)
  coldata$time <- as.factor(coldata$time)
  all(rownames(coldata) %in% colnames(cts2))
  cts2 <- cts2[, rownames(coldata)]
  all(rownames(coldata) == colnames(cts2))
  data <- list(cts2, coldata)
  return(data)  
}

get_deseq <- function(cts, coldata, sample_list, rowsums, thr){
  dds <- DESeqDataSetFromMatrix(countData = cts[,sample_list],
                              colData = coldata[sample_list,],
                              design = ~ condition)
                              #design = ~ time)
  print(nrow(dds))
  dds<-dds[rowSums( counts(dds) != 0 ) >= rowsums,] 
  dds<-dds[rowSums(counts(dds)) >= thr,] 
  print(nrow(dds))
  dds$condition <- relevel(dds$condition, ref = "C")
  dds <- DESeq(dds)
  print(nrow(dds))
  return(dds)
}

get_deseq_ref_D <- function(cts, coldata, sample_list, rowsums){
  dds <- DESeqDataSetFromMatrix(countData = cts[,sample_list],
                                colData = coldata[sample_list,],
                                design = ~ condition)
  print(nrow(dds))
  dds<-dds[rowSums( counts(dds) != 0 ) >= rowsums,] 
  print(nrow(dds))
  dds$condition <- relevel(dds$condition, ref = "D")
  dds <- DESeq(dds)
  print(nrow(dds))
  return(dds)
}

get_deseq_mp <- function(cts, coldata, sample_list, rowsums){
  dds <- DESeqDataSetFromMatrix(countData = cts[,sample_list],
                                colData = coldata[sample_list,],
                                design = ~ batch + condition)
  #design = ~ time)
  print(nrow(dds))
  dds<-dds[rowSums( counts(dds) != 0 ) >= rowsums,] 
  print(nrow(dds))
  dds$condition <- relevel(dds$condition, ref = "C")
  dds <- DESeq(dds)
  print(nrow(dds))
  return(dds)
}

get_deseq_batch <- function(cts, coldata, sample_list, rowsums,thr){
  dds <- DESeqDataSetFromMatrix(countData = cts[,sample_list],
                                colData = coldata[sample_list,],
                                design = ~ block + condition)
  #design = ~ time)
  print(nrow(dds))
  dds<-dds[rowSums( counts(dds) != 0 ) >= rowsums,] 
  print(nrow(dds))
  dds<-dds[rowSums(counts(dds)) >= thr,] 
  print(nrow(dds))
  dds$condition <- relevel(dds$condition, ref = "C")
  dds <- DESeq(dds, test="LRT", reduced=~block)
  #dds <- DESeq(dds)
  print(nrow(dds))
  return(dds)
}

get_deseq_blood_time <- function(cts, coldata, sample_list, rowsums){
  dds <- DESeqDataSetFromMatrix(countData = cts[,sample_list],
                                colData = coldata[sample_list,],
                                design = ~ subject + time)
  print(nrow(dds))
  dds<-dds[rowSums( counts(dds) != 0 ) >= rowsums,] 
  print(nrow(dds))
  dds$time <- relevel(dds$time, ref = "0")
  dds <- DESeq(dds, test="LRT", reduced=~subject)
  #dds <- DESeq(dds)
  return(dds)
}

# draw basic data exploration plots - heatmap, sample distance, PCA. 
# vsd can be replaced with other normalization results - e.g. rld, ntd   
basic_plots <-function(dds, vsd, org, tissue, type, time=FALSE){ 
  #Mean - SD relationship
  png(paste("figures/meanSD_", org, "_", tissue, ".png", sep=""))
  meanSdPlot(assay(vsd))
  dev.off()
  #Dispersion estimate vs mean
  png(paste("figures/dispEst_", org, "_", tissue, ".png", sep=""))
  plotDispEsts(dds)
  dev.off()
  #Cook distance boxplot
  png(paste("figures/CooksDist_", org, "_", tissue, ".png", sep=""))
  par(mar=c(8,5,2,2))
  boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
  dev.off()
  #Heatmap
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:20]
  df <- as.data.frame(colData(dds)[,c("condition","time")])
  #df <- as.data.frame(colData(dds)[,c("condition","subject")])
  #print(df)
  x <- pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
              cluster_cols=TRUE, annotation_col=df)
  save_pheatmap_pdf(x, paste("figures/heatmap_", org, "_",tissue,"_",type,".pdf", sep="")) 
  #Samples-distance plot
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rownames(colData(dds)), vsd$condition, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  x<-pheatmap(sampleDistMatrix,
              clustering_distance_rows=sampleDists,
              clustering_distance_cols=sampleDists,
              fontsize_row = 20,
              col=colors)
  save_pheatmap_pdf(x, paste("figures/sample_distance_", org, "_",tissue,"_",type,".pdf", sep=""),width=9,height=7) 
  if (time==TRUE) {group="time"}
  else {group="condition"}
  print(group)
  plotPCA(vsd, intgroup=c(group)) +
    #geom_point(aes(size=25))+
    aes(label = name) +
    geom_text_repel()+
    theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"),
          plot.title = element_text(hjust = 0.5,size=20,face="bold"),legend.text = element_text(size=12, face="bold"))
  ggsave(paste("figures/PCA_", org, "_",tissue,"_",type, ".png", sep=""), width=5, height=3.5, dpi = 1000)
  
}

get_results <- function(dds, org, tissue, type, thr, orth1="mm", orth2="pp",biotypes,timepoint){
  if (tissue=="colon"){
    res <- results(dds, independentFiltering = FALSE, contrast=c("condition","D","C"))
  }
  if (tissue=="blood"){
    res <- results(dds, independentFiltering = FALSE, contrast=c("time",timepoint,"0"))
    }
  resOrdered <- as.data.frame(res[order(abs(res$log2FoldChange)),])
  names(resOrdered)[names(resOrdered) == "log2FoldChange"] = "logFC"
  names(resOrdered)[names(resOrdered) == "baseMean"] = "avg.expr"
  if (type == "small"){ 
    resOrdered$id <-tolower(gsub("_.*","",rownames(resOrdered)))
  }
  else if (type == "total"){
    resOrdered$id <- rownames(resOrdered)
    resOrdered$id <- gsub("\\..*","",resOrdered$id)
    resOrdered <- as.data.frame(merge(resOrdered, orth1, by.x='id', by.y='Gene.stable.ID', all.x=TRUE))
    resOrdered <- as.data.frame(merge(resOrdered, orth2, by.x='id', by.y='Gene.stable.ID', all.x=TRUE))
  }
  write.table(resOrdered[c(1,2,6,7,8,9)], file=paste("results/", org, "_", tissue, "_",type, "_", thr, "_thr", ".tsv", sep=""), sep="\t")
  #resOrdered_uc_co<-resOrdered_uc_co[resOrdered_uc_co$padj < thr,]
  resOrdered <- merge(resOrdered, biotypes, by.x='id',by.y='id')
  resSig <- subset(resOrdered, padj <= thr)
  return(resSig)
}

pca_plot <- function(vsd, org, tissue, type){
  #org="mouse"; tissue="colon"; type="small"
  plotPCA(vsd, intgroup=c("condition")) + 
    aes(label = name,size=9)+
    geom_text_repel()+
    theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),plot.title = element_text(hjust = 0.5,size=20,face="bold"))
  ggsave(paste("figures/PCA_", org, "_",tissue,"_",type, ".png", sep=""), width=9, height=7, dpi = 100)

}

plot_gene <- function(gene,gene_name,dds){
  #https://github.com/hbctraining/DGE_workshop/blob/master/lessons/06_DGE_visualizing_results.md
  #pdf(file=paste("figures/",gene,".pdf",sep="")) 
  p <- plotCounts(dds, gene=gene, intgroup="condition", returnData=TRUE)
  # Plotting the gene normalized counts, using the samplenames (rownames of d as labels)
  p1<-ggplot(p, aes(x = condition, y = count, color = condition, fill = condition)) + 
    geom_dotplot(binaxis="y", stackdir="center")+
    aes(label = rownames(p)) +
    geom_text_repel() +
    ggtitle(paste(gene,gene_name,sep=", ")) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste("figures/", gene_name, ".png", sep=""), width=5, height=4, dpi = 1000)
  plot(p1)
} 

plot_two_genes <- function(gene1, gene_name1, dds1, gene2, gene_name2, dds2){
  #https://github.com/hbctraining/DGE_workshop/blob/master/lessons/06_DGE_visualizing_results.md
  #pdf(file=paste("figures/",gene,".pdf",sep="")) 
  if (gene1!=""){
  p1 <- plotCounts(dds1, gene=gene1, intgroup="condition", returnData=TRUE)
  # Plotting the gene normalized counts, using the samplenames (rownames of d as labels)
  p11 <- ggplot(p1, aes(x = condition, y = count, fill = condition, color=condition)) + 
    #geom_point(position=position_jitter(w = 0.05,h = 0)) +
    geom_dotplot(binaxis="y", stackdir="center")+
    expand_limits(y=0) +
    aes(label = rownames(p1)) +
    geom_text_repel() +
    ggtitle(paste(gene1,gene_name1,sep=", ")) +
    theme(plot.title = element_text(hjust = 0.5))
  }
  if (gene2!=""){
  p2 <- plotCounts(dds2, gene=gene2, intgroup="condition", returnData=TRUE)
  # Plotting the gene normalized counts, using the samplenames (rownames of d as labels)
  p22 <- ggplot(p2, aes(x = condition, y = count, color = condition, fill=condition)) + 
    geom_dotplot(binaxis="y", stackdir="center")+
    expand_limits(y=0) +
    aes(label = rownames(p2)) +
    geom_text_repel() +
    ggtitle(paste(gene2,gene_name2,sep=", ")) +
    theme(plot.title = element_text(hjust = 0.5))
  }
  if (gene2!="" & gene1!=""){ 
    p<-cowplot::plot_grid(p11, p22, labels = "AUTO")
  }
  else if (gene1==""){
    p<-cowplot::plot_grid(NULL, p22, labels = "AUTO")
  }
  else if (gene2==""){
    p<-cowplot::plot_grid(p11, NULL, labels = "AUTO")
  }
  plot(p)
  #ggsave(paste("figures/", gene1,"_", gene2, ".png", sep=""), width=8, height=3, dpi = 1000)
  ggsave(paste("figures/", gene_name1, ".png", sep=""),p, width=8, height=3, dpi = 1000)
  #plot(p)
} 

plot_three_genes <- function(gene1, gene_name1, dds1, gene2, gene_name2, dds2, array_expr, groups){
  #https://github.com/hbctraining/DGE_workshop/blob/master/lessons/06_DGE_visualizing_results.md
  #pdf(file=paste("figures/",gene,".pdf",sep="")) 
  if (gene1!=""){
    p1 <- plotCounts(dds1, gene=gene1, intgroup="condition", returnData=TRUE)
    # Plotting the gene normalized counts, using the samplenames (rownames of d as labels)
    p11 <- ggplot(p1, aes(x = condition, y = count, fill = condition, color=condition)) + 
      #geom_point(position=position_jitter(w = 0.05,h = 0)) +
      geom_dotplot(binaxis="y", stackdir="center")+
      aes(label = rownames(p1)) +
      geom_text_repel() +
      ggtitle(paste(gene1,gene_name1,sep=", ")) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  if (gene2!=""){
    p2 <- plotCounts(dds2, gene=gene2, intgroup="condition", returnData=TRUE)
    # Plotting the gene normalized counts, using the samplenames (rownames of d as labels)
    p22 <- ggplot(p2, aes(x = condition, y = count, color = condition, fill=condition)) + 
      geom_dotplot(binaxis="y", stackdir="center")+
      aes(label = rownames(p2)) +
      geom_text_repel() +
      ggtitle(paste(gene2,gene_name2,sep=", ")) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  a<-plot_array_int(gene_name1, array_expr, groups)
  #a<-plot(jitter(as.integer(groups)),array_expr, ylab="log intensity", xaxt="n",
  #        xlab="",xlim=c(0.5, nlevels(groups)+0.5), main=gene_name1)
  #a+axis(side=1, at=seq_len(nlevels(groups)), label=levels(groups))
  if (gene2!="" & gene1!=""){ 
    #p<-cowplot::plot_grid(p11, p22, a, labels = "AUTO")
    p<-cowplot::plot_grid(p11,  p22, labels = "AUTO")
  }
  else if (gene1==""){
    p<-cowplot::plot_grid(NULL, p22, a, labels = "AUTO")
  }
  else if (gene2==""){
    p<-cowplot::plot_grid(p11, a , NULL, labels = "AUTO")
  }
  plot(p)
  #ggsave(paste("figures/", gene1,"_", gene2, ".png", sep=""), width=8, height=3, dpi = 1000)
  ggsave(paste("figures/", gene_name1, ".png", sep=""),p, width=8, height=3, dpi = 1000)
  #plot(p)
  plot(a)
} 

plot_array_int <- function(gene_name1, array_expr, groups){
  a<-plot(jitter(as.integer(groups)),array_expr, ylab="log intensity", xaxt="n",
          xlab="", xlim=c(0.5, nlevels(groups)+0.5), main=gene_name1)
  a+axis(side=1, at=seq_len(nlevels(groups)), label=levels(groups))
  return(a)
}

get_orth_table <- function(orth_df, test_col, keep_cols){
  orth_o2o <- orth_df[orth_df[,test_col]=="ortholog_one2one",][,keep_cols]
  orth_o2o <- orth_o2o[!duplicated(orth_o2o), ]
  return(orth_o2o)
}

get_orth_table_all <- function(orth_df, test_col, keep_cols){
  orth_o2o <- orth_df[orth_df[,test_col]=="ortholog_one2one",][,keep_cols]
  orth_o2o <- orth_o2o[!duplicated(orth_o2o), ]
  return(orth_o2o)
}

get_correlations_pct <- function(common, species1, species2) {
  ap = cor(common$qAI, common$qCD)
  bp = cor(common$qAI, common$qCC)
  cp = cor(common$qC, common$qCC)
  dp = cor(common$qC, common$qCD)
  as = cor(common$qAI, common$qCD, method="spearman")
  bs = cor(common$qAI, common$qCC, method="spearman")
  cs = cor(common$qC, common$qCC, method="spearman")
  ds = cor(common$qC, common$qCD, method="spearman")
  pears <- matrix(c(ap,dp,bp,cp), nrow = 2, dimnames = list(c(paste(species1,"_disease",sep=""),paste(species1,"_healthy",sep="")), c(paste(species2,"_disease",sep=""),paste(species2,"_healthy",sep=""))))
  spear <- matrix(c(as,ds,bs,cs), nrow = 2, dimnames = list(c(paste(species1,"_disease",sep=""),paste(species1,"_healthy",sep="")), c(paste(species2,"_disease",sep=""),paste(species2,"_healthy",sep=""))))
  return(pears)
  #return(spear)
}

get_correlations <- function(common, species1, species2) {
  ap = cor(common$AI, common$CD)
  bp = cor(common$AI, common$CC)
  cp = cor(common$C, common$CC)
  dp = cor(common$C, common$CD)
  as = cor(common$AI, common$CD, method="spearman")
  bs = cor(common$AI, common$CC, method="spearman")
  cs = cor(common$C, common$CC, method="spearman")
  ds = cor(common$C, common$CD, method="spearman")
  pears <- matrix(c(ap,dp,bp,cp), nrow = 2, dimnames = list(c(paste(species1,"_disease",sep=""),paste(species1,"_healthy",sep="")), c(paste(species2,"_disease",sep=""),paste(species2,"_healthy",sep=""))))
  spear <- matrix(c(as,ds,bs,cs), nrow = 2, dimnames = list(c(paste(species1,"_disease",sep=""),paste(species1,"_healthy",sep="")), c(paste(species2,"_disease",sep=""),paste(species2,"_healthy",sep=""))))
  return(pears)
  #return(spear)
}

get_correlations_blood <- function(common, species1, species2) {
  ap = cor(common$AI, common$CD)
  bp = cor(common$AI, common$CC)
  cp = cor(common$C, common$CC)
  dp = cor(common$C, common$CD)
  as = cor(common$AI, common$CD, method="spearman")
  bs = cor(common$AI, common$CC, method="spearman")
  cs = cor(common$C, common$CC, method="spearman")
  ds = cor(common$C, common$CD, method="spearman")
  pears <- matrix(c(ap,dp,bp,cp), nrow = 2, dimnames = list(c(paste(species1,"_disease",sep=""),paste(species1,"_healthy",sep="")), c(paste(species2,"_disease",sep=""),paste(species2,"_healthy",sep=""))))
  spear <- matrix(c(as,ds,bs,cs), nrow = 2, dimnames = list(c(paste(species1,"_disease",sep=""),paste(species1,"_healthy",sep="")), c(paste(species2,"_disease",sep=""),paste(species2,"_healthy",sep=""))))
  return(pears)
  #return(spear)
}

addmean <- function(df, samples, name){
  df<-cbind(df,rowMeans(df[,samples]))
  colnames(df)[ncol(df)] <- name
  return(df)
}

transform_tpm <- function(tpm, coldata, species){
  tpm<-tpm[rowSums( tpm != 0 ) >= 3,] 
  tpm <- log(tpm+1)
  CC <- rownames(coldata[coldata$tissue=="colon" & coldata$condition=="C",]) 
  CD <- rownames(coldata[coldata$tissue=="colon" & coldata$condition=="D",]) 
  CA <- rownames(coldata[coldata$tissue=="colon",])
  BC <- rownames(coldata[coldata$tissue=="blood" & coldata$time=="1" & coldata$condition=="D",])
  BC <- rownames(coldata[coldata$tissue=="blood" & coldata$condition=="C",])
  #BC <- rownames(coldata[coldata$tissue=="blood" & coldata$time=="1",])
  if (species=="mouse"){
    BD <- rownames(coldata[coldata$tissue=="blood" & coldata$condition=="D" & coldata$time=="3",])}
  if (species=="pig"){
    BD <- rownames(coldata[coldata$tissue=="blood" & coldata$condition=="D" & coldata$time=="4",])}
  BA <- rownames(coldata[coldata$tissue=="blood",])
  tpm <- addmean(tpm, CC, "CC")
  tpm <- addmean(tpm, CD, "CD")
  tpm <- addmean(tpm, CA, "CA")
  tpm <- addmean(tpm, BC, "BC")
  tpm <- addmean(tpm, BD, "BD")
  tpm <- addmean(tpm, BA, "BA")
  tpm <- as.data.frame(tpm)
  #tpm <- within(tpm, qCC <- .bincode(CC, quantile(CA, probs=0:100/100), include.lowest=TRUE))
  #tpm <- within(tpm, qCD <- .bincode(CD, quantile(CA, probs=0:100/100), include.lowest=TRUE))
  #tpm <- within(tpm, qBC <- .bincode(BC, quantile(BA, probs=0:100/100), include.lowest=TRUE))
  #tpm <- within(tpm, qBD <- .bincode(BD, quantile(BA, probs=0:100/100), include.lowest=TRUE))
  tpm$id <- rownames(tpm)
  if (species=="mouse"){
    tpm <- as.data.frame(merge(tpm, mh_orth, by.x='id', by.y='Gene.stable.ID'))#, all.x=TRUE))
    #tpm <- as.data.frame(merge(tpm, mp_orth, by.x='id', by.y='Gene.stable.ID', all.x=TRUE))
  }
  if (species=="pig"){
    tpm <- as.data.frame(merge(tpm, ph_orth, by.x='id', by.y='Gene.stable.ID'))#, all.x=TRUE))
    #tpm <- as.data.frame(merge(tpm, pm_orth, by.x='id', by.y='Gene.stable.ID', all.x=TRUE))
  }
  return(tpm)
}

draw_logFC_corr <- function(data1, data2, common, p1, p2){
  pdf(file=paste("figures/logFC_correlation_",data1,"_",data2,".pdf",sep=""),width=6,height=5) 
  plot(common$logFC.x, common$logFC.y, col="darkgrey", panel.first=grid(),
       main=paste(data1,"_",data2, sep=''), xlab=paste("logFC ",data1,sep=''), ylab=paste("logFC ",data2, sep=''),
       pch=20, cex=0.6)
  abline(v=0)
  abline(h=0)
  with(subset(common, padj.x<p1 ), points(logFC.x, logFC.y, pch=20, col="red", cex=0.7))
  with(subset(common, padj.y<p2 ), points(logFC.x, logFC.y, cex=1.1))
  selected <- subset(common, abs(common$logFC.x)>1.5 & abs(common$logFC.y)>1.5) #& common$logFC.x*common$logFC.y > 0
  #print(selected)
  text(selected$logFC.x, selected$logFC.y,pos=3,
       lab=selected$Gene.Symbol, cex=0.4)
  dev.off()
}

get_stats <- function(common, species1, species2, pthr1,pthr2, lfcthr1, lfcthr2, exprthr){
  sig1 <- common[abs(common$logFC.x)>lfcthr1 & common$padj.x<pthr1, ]
  sig2 <- common[abs(common$logFC.y)>lfcthr2 & common$padj.y<pthr2, ]
  sig <- common[abs(common$logFC.x)>lfcthr1 & abs(common$logFC.y)>lfcthr2 & 
                  common$padj.x<pthr1 & common$padj.y<pthr2 & 
                  common$avg.expr.y > exprthr &
                  common$logFC.x*common$logFC.y > 0, ]
  opposite <- common[abs(common$logFC.x)>lfcthr1 & abs(common$logFC.y)>lfcthr2 & 
                       common$padj.x<pthr1 & common$padj.y<pthr2 & 
                       common$logFC.x*common$logFC.y < 0, ]
  stats <- c(dim(common)[1], dim(sig1)[1], dim(sig2)[1], dim(sig)[1], dim(opposite)[1])
  #dimnames = list(c(species1, species2, "consistent dysreg", "inconsistent dysreg"),c(""))
  out <- matrix(stats, nrow = 1)
  write.table(sig, file=paste("results/consistent_",species1,"_",species2, ".tsv",sep=""), row.names=F, sep="\t")
  write.table(opposite, file=paste("results/inconsistent_",species1,"_",species2, ".tsv",sep=""), row.names=F, sep="\t")
  colnames(out) <- c("total", species1, species2, "consistent", "opposite" )
  data <- list(out, sig, opposite)
  return(data)
}
