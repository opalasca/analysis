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
library(knitr) 
library("fdrtool") 

#require("twbattaglia/btools")

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
  coldata$condition <- as.character(coldata$condition)
  coldata$condition[coldata$condition == "D"] <- "DSS"
  coldata$condition[coldata$condition == "C"] <- "Control"
  coldata$condition <- as.factor(coldata$condition)
  if ("block" %in% colnames(coldata)){
    coldata$block <- as.factor(coldata$block)
  }
  print(all(rownames(coldata) %in% colnames(cts)))
  cts <- cts[, rownames(coldata)]
  print(all(rownames(coldata) == colnames(cts)))
  data <- list(cts, coldata)
  return(data)  
}

get_data_miRNA_v2 <- function(counts_file, config_file, n, selected){
  cts <- data.frame()
  coldata <- data.frame()
  cts <- read.csv(counts_file, check.names=FALSE, header=TRUE, sep="\t")
  colnames(cts)[1]<-"miRNA"
  #cts$miRNA <-tolower(gsub("_.*","",cts$miRNA))
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
  coldata$condition <- as.character(coldata$condition)
  coldata$condition[coldata$condition == "D"] <- "DSS"
  coldata$condition[coldata$condition == "C"] <- "Control"
  coldata$condition <- as.factor(coldata$condition)
  if ("block" %in% colnames(coldata)){
    coldata$block <- as.factor(coldata$block)
  }
  
  all(rownames(coldata) %in% colnames(cts2))
  cts2 <- cts2[, rownames(coldata)]
  all(rownames(coldata) == colnames(cts2))
  print(nrow(cts2))
  cts2 <- cts2[rownames(cts2) %in% selected,]
  print(nrow(cts2))
  data <- list(cts2, coldata)
  return(data)  
}

get_data_miRNA <- function(counts_file, config_file, n){
  cts <- data.frame()
  coldata <- data.frame()
  cts <- read.csv(counts_file, check.names=FALSE, header=TRUE, sep="\t")
  colnames(cts)[1]<-"miRNA"
  cts$miRNA <-tolower(gsub("_.*","",cts$miRNA))
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
  coldata$condition <- as.character(coldata$condition)
  coldata$condition[coldata$condition == "D"] <- "DSS"
  coldata$condition[coldata$condition == "C"] <- "Control"
  coldata$condition <- as.factor(coldata$condition)
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
  dds <- estimateSizeFactors(dds)
  mnc <- rowMeans(counts(dds, normalized=TRUE))
  dds <- dds[mnc > 5,]
  #dds<-dds[rowSums( counts(dds) != 0 ) >= rowsums,] 
  #dds<-dds[rowSums(counts(dds)) >= thr,] 
  print(nrow(dds))
  dds$condition <- relevel(dds$condition, ref = "Control")
  dds <- DESeq(dds)
  print(nrow(dds))
  return(dds)
}

get_deseq_time_cond_merged <- function(cts, coldata, sample_list, rowsums, thr, refcond){
  # A design that takes all blood samples into consideration, by treating all control timepoints,
  # together with day0DSS as control
  coldata$timecond <- paste("day",coldata$time,"_",coldata$condition,sep="")
  coldata$timecond[coldata$condition=="Control"] <- "day0_DSS" 
  coldata$timecond <- as.factor(coldata$timecond)
  dds <- DESeqDataSetFromMatrix(countData = cts[,sample_list],
                                colData = coldata[sample_list,],
                                design = ~ timecond)
  print(nrow(dds))
  dds <- estimateSizeFactors(dds)
  mnc <- rowMeans(counts(dds, normalized=TRUE))
  dds <- dds[mnc > 3,]
  print(nrow(dds))
  dds$timecond <- relevel(dds$timecond, ref = refcond)
  dds <- DESeq(dds)
  print(nrow(dds))
  return(dds)
}

get_deseq_time_cond_merged_plus_batch <- function(cts, coldata, sample_list, rowsums, thr, refcond){
  # A design that takes all blood samples into consideration, by treating all control timepoints,
  # together with day0DSS as control
  # subject batch effect needed for pig - helps a lot with the p val distribution
  coldata$timecond <- paste("day",coldata$time,"_",coldata$condition,sep="")
  coldata$timecond[coldata$condition=="Control"] <- "day0_DSS" 
  coldata$timecond <- as.factor(coldata$timecond)
  dds <- DESeqDataSetFromMatrix(countData = cts[,sample_list],
                                colData = coldata[sample_list,],
                                design = ~ subject + timecond)
  print(nrow(dds))
  dds <- estimateSizeFactors(dds)
  mnc <- rowMeans(counts(dds, normalized=TRUE))
  dds <- dds[mnc > 3,]
  print(nrow(dds))
  dds$timecond <- relevel(dds$timecond, ref = refcond)
  print(colData(dds))
  dds <- DESeq(dds, test="LRT", reduced=~subject)
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
  dds$condition <- relevel(dds$condition, ref = "DSS")
  dds <- DESeq(dds)
  print(nrow(dds))
  return(dds)
}

get_deseq_batch <- function(cts, coldata, sample_list, rowsums,thr){
  dds <- DESeqDataSetFromMatrix(countData = cts[,sample_list],
                                colData = coldata[sample_list,],
                                design = ~ block + condition)
  print(nrow(dds))
  dds<-dds[rowSums( counts(dds) != 0 ) >= rowsums,] 
  dds <- estimateSizeFactors(dds)
  mnc <- rowMeans(counts(dds, normalized=TRUE))
  dds <- dds[mnc > 5,]
  print(nrow(dds))
  dds<-dds[rowSums(counts(dds)) >= thr,] 
  print(nrow(dds))
  dds$condition <- relevel(dds$condition, ref = "Control")
  dds <- DESeq(dds, test="LRT", reduced=~block)
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
  return(dds)
}

get_deseq_blood_time_no_subject <- function(cts, coldata, sample_list, rowsums){
  dds <- DESeqDataSetFromMatrix(countData = cts[,sample_list],
                                colData = coldata[sample_list,],
                                design = ~ time)
  print(nrow(dds))
  dds<-dds[rowSums( counts(dds) != 0 ) >= rowsums,] 
  print(nrow(dds))
  dds$time <- relevel(dds$time, ref = "0")
  dds <- DESeq(dds)
  return(dds)
}

get_deseq_blood_time_full <- function(cts, coldata, sample_list, rowsums){
  dds <- DESeqDataSetFromMatrix(countData = cts[,sample_list],
                                colData = coldata[sample_list,],
                                design = ~ condition + time + condition:time)
  print(nrow(dds))
  dds<-dds[rowSums( counts(dds) != 0 ) >= rowsums,] 
  print(nrow(dds))
  dds$time <- relevel(dds$time, ref = "0")
  dds <- DESeq(dds, test="LRT", reduced = ~ condition + time)
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
  #save_pheatmap_pdf(x, paste("figures/heatmap_", org, "_",tissue,"_",type,".pdf", sep="")) 
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
  pdf("figures/PCA_3D.pdf")
  plotPCA3D(vsd, intgroup = "condition", ntop = 500,
            returnData = FALSE)
  dev.off()
  #ggsave(paste("figures/PCA_3D_", org, "_",tissue,"_",type, ".png", sep=""), width=9, height=7, dpi = 100)
}

get_results <- function(dds, org, tissue, type, thr, orth1="mm", orth2="pp", biotypes, idname, mir_seq=NULL, timepoint, cond1, cond2){
  if (timepoint=="" & is.null(cond1)){
    print("no timepoint")
    res <- results(dds,cooksCutoff=FALSE, independentFiltering = FALSE, test="Wald",contrast=c("condition","DSS","Control"))
  }
  if (tissue=="blood" & timepoint!="" & is.null(cond1)){
    print("time given")
    res <- results(dds, cooksCutoff=FALSE, independentFiltering = FALSE, test="Wald", contrast=c("time",timepoint,"0"))
  }
  if (!is.null(cond1) & !is.null(cond2)){
    print("time and condition merged")
    res <- results(dds, cooksCutoff=FALSE, independentFiltering = FALSE, test="Wald", contrast=c("timecond",cond1,cond2))
  }
  resOrdered <- as.data.frame(res[order(abs(res$log2FoldChange)),])
  names(resOrdered)[names(resOrdered) == "log2FoldChange"] = "logFC"
  names(resOrdered)[names(resOrdered) == "baseMean"] = "avg.expr"
  if (type == "small"){ 
    resOrdered$id <- rownames(resOrdered)
    resOrdered$partial_id <- gsub("^[^-]*-","",resOrdered$id)
    resOrdered <- as.data.frame(resOrdered)
    resOrdered$biotype <- "miRNA"
    if (!is.null(mir_seq)){
      resOrdered <- as.data.frame(merge(resOrdered, mir_seq, by.x='id', by.y='id', all.x=TRUE))
    }
    write.table(resOrdered[c(1,2,3,6,7,8,10)], file=paste("results/", org, "_", tissue, "_",type, "_", thr, "_thr", ".tsv", sep=""), sep="\t", row.names = FALSE, quote = F)
  }
  else if (startsWith(type, "total")){
    idname <- idname[,c(1,3)]
    names(idname) <- c("id","Gene.Name")
    idname <- unique(idname)
    resOrdered$id <- rownames(resOrdered)
    resOrdered$id <- gsub("\\..*","",resOrdered$id)
    resOrdered <- as.data.frame(merge(resOrdered, orth1, by.x='id', by.y='Gene.stable.ID', all.x=TRUE))
    resOrdered <- as.data.frame(merge(resOrdered, orth2, by.x='id', by.y='Gene.stable.ID', all.x=TRUE))
    resOrdered <- as.data.frame(merge(resOrdered, idname, by.x='id', by.y='id', all.x=TRUE))
    write.table(resOrdered[c(1,2,3,6,7,12,8,9,10,11)], file=paste("results/", org, "_", tissue, "_",type, "_", thr, "_thr", ".tsv", sep=""), sep="\t", row.names = FALSE, quote = F)
    resOrdered <- merge(resOrdered, biotypes, by.x='id',by.y='id', all.x=TRUE)
}
  #resOrdered_uc_co<-resOrdered_uc_co[resOrdered_uc_co$padj < thr,]
  resSig <- subset(resOrdered, padj <= thr)
  return(resSig)
}

get_results_corrected_pval <- function(dds, org, tissue, type, thr, orth1="mm", orth2="pp", biotypes, mir_seq=NULL, timepoint, cond1, cond2){
  if (timepoint=="" & is.null(cond1)){
    print("no timepoint")
    cooksCutoff=FALSE
    res <- results(dds, cooksCutoff=FALSE, independentFiltering = FALSE, test="Wald",contrast=c("condition","DSS","Control"))
  }
  if (tissue=="blood" & timepoint!="" & is.null(cond1)){
    print("time given")
    res <- results(dds, cooksCutoff=FALSE, independentFiltering = FALSE, test="Wald", contrast=c("time",timepoint,"0"))
  }
  if (!is.null(cond1) & !is.null(cond2)){
    print("time and condition merged")
    res <- results(dds, cooksCutoff=FALSE, independentFiltering = FALSE, test="Wald", contrast=c("timecond",cond1,cond2))
  }
  
  DESeq2Res<-res
  #DESeq2Res <- DESeq2Res[, -which(names(DESeq2Res) == "padj")]
  FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = F)
  FDR.DESeq2Res$param[1, "sd"]
  DESeq2Res[,"pvalue"]  <- FDR.DESeq2Res$pval
  DESeq2Res[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
  res<- as.data.frame(DESeq2Res)
  
  resOrdered <- as.data.frame(res[order(abs(res$log2FoldChange)),])
  names(resOrdered)[names(resOrdered) == "log2FoldChange"] = "logFC"
  names(resOrdered)[names(resOrdered) == "baseMean"] = "avg.expr"
  if (type == "small"){ 
    resOrdered$id <- rownames(resOrdered)
    resOrdered$partial_id <- gsub("^[^-]*-","",resOrdered$id)
    resOrdered <- as.data.frame(resOrdered)
    resOrdered$biotype <- "miRNA"
    if (!is.null(mir_seq)){
      resOrdered <- as.data.frame(merge(resOrdered, mir_seq, by.x='id', by.y='id', all.x=TRUE))
    }
    write.table(resOrdered[c(1,2,3,6,7,8,10)], file=paste("results/", org, "_", tissue, "_",type, "_", thr, "_thr", ".tsv", sep=""), sep="\t", row.names = FALSE, quote = F)
  }
  else if (startsWith(type, "total")){
    resOrdered$id <- rownames(resOrdered)
    resOrdered$id <- gsub("\\..*","",resOrdered$id)
    resOrdered <- as.data.frame(merge(resOrdered, orth1, by.x='id', by.y='Gene.stable.ID', all.x=TRUE))
    resOrdered <- as.data.frame(merge(resOrdered, orth2, by.x='id', by.y='Gene.stable.ID', all.x=TRUE))
    write.table(resOrdered[c(1,2,3,7,8,9)], file=paste("results/", org, "_", tissue, "_",type, "_", thr, "_thr", ".tsv", sep=""), sep="\t", row.names = FALSE, quote = F)
    resOrdered <- merge(resOrdered, biotypes, by.x='id',by.y='id', all.x=TRUE)
  }
  #resOrdered_uc_co<-resOrdered_uc_co[resOrdered_uc_co$padj < thr,]
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
  plotPCA3D(vsd, intgroup = "condition", ntop = 500,
            returnData = FALSE)
  ggsave(paste("figures/PCA_3D_", org, "_",tissue,"_",type, ".png", sep=""), width=9, height=7, dpi = 100)
  
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
  CC <- rownames(coldata[coldata$tissue=="colon" & coldata$condition=="Control",]) 
  CD <- rownames(coldata[coldata$tissue=="colon" & coldata$condition=="DSS",]) 
  CA <- rownames(coldata[coldata$tissue=="colon",])
  BC <- rownames(coldata[coldata$tissue=="blood" & coldata$time=="1" & coldata$condition=="DSS",])
  BC <- rownames(coldata[coldata$tissue=="blood" & coldata$condition=="Control",])
  #BC <- rownames(coldata[coldata$tissue=="blood" & coldata$time=="1",])
  if (species=="mouse"){
    BD <- rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" & coldata$time=="3",])}
  if (species=="pig"){
    BD <- rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" & coldata$time=="4",])}
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

get_stats_pval <- function(common, species1, species2, pthr1,pthr2, lfcthr1, lfcthr2, exprthr){
  sig1 <- common[abs(common$logFC.x)>lfcthr1 & common$pvalue.x<pthr1, ]
  sig2 <- common[abs(common$logFC.y)>lfcthr2 & common$pvalue.y<pthr2, ]
  sig <- common[abs(common$logFC.x)>lfcthr1 & abs(common$logFC.y)>lfcthr2 & 
                  common$pvalue.x<pthr1 & common$pvalue.y<pthr2 & 
                  common$avg.expr.y > exprthr &
                  common$logFC.x*common$logFC.y > 0, ]
  opposite <- common[abs(common$logFC.x)>lfcthr1 & abs(common$logFC.y)>lfcthr2 & 
                       common$pvalue.x<pthr1 & common$pvalue.y<pthr2 & 
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

heatmap_DE <- function(sig, rld, org, tissue, type, title, timepoint){
  rows <- match(sig$id, row.names(rld))
  mat <- assay(rld)[rows,]
  n=dim(mat)[[1]]
  df <- as.data.frame(colData(rld))
  df = df[-c(2:7)]
  #breaksList = seq(-2, 2, by = 0.25)
  annot_heatmap = c("DSS1","DSS2","DSS3","Control1","Control2","Control3")
  x <- pheatmap(mat, 
                scale="row",
                cluster_rows=TRUE,
                cluster_cols =FALSE, 
                show_rownames=FALSE,
                #show_colnames=FALSE,
                annotation_col=df,
                annotation_colors = list(condition=c(DSS="grey30", Control="grey70")),
                #annotation_legend = FALSE,
                annotation_names_col = FALSE,
                #annotation_col=c("DSS","Control"),
                key.title = "Row-wise Z-score",
                #breaks = breaksList,
                border_color = "NA",
                #color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
                #cellwidth=20,
                main=paste(title,",n=",n,sep=""))
  save_pheatmap_pdf(x, paste("figures/heatmap_DE_", org, "_", tissue, "_", type, "_", timepoint, ".pdf", sep=""))
  #,width=9,height=7
  clust <- cbind(mat, cluster = cutree(x$tree_row, k = 5))
  return(clust)
}

heatmap_DE_rownames <- function(sig, rld, org, tissue, type, title, timepoint){
  rows <- match(sig$id, row.names(rld))
  mat <- assay(rld)[rows,]
  n=dim(mat)[[1]]
  df <- as.data.frame(colData(rld))
  df = df[-c(2:7)]
  #breaksList = seq(-2, 2, by = 0.25)
  annot_heatmap = c("DSS1","DSS2","DSS3","Control1","Control2","Control3")
  x <- pheatmap(mat, 
                scale="row",
                cluster_rows=TRUE,
                cluster_cols =FALSE, 
                #show_rownames=FALSE,
                #show_colnames=FALSE,
                annotation_col=df,
                annotation_colors = list(condition=c(DSS="grey30", Control="grey70")),
                #annotation_legend = FALSE,
                annotation_names_col = FALSE,
                #annotation_col=c("DSS","Control"),
                key.title = "Row-wise Z-score",
                #breaks = breaksList,
                border_color = "NA",
                #color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
                #cellwidth=20,
                main=paste(title,",n=",n,sep=""))
  save_pheatmap_pdf(x, paste("figures/heatmap_DE_", org, "_", tissue, "_", type, "_", timepoint, ".pdf", sep=""))
  #,width=9,height=7
  clust <- cbind(mat, cluster = cutree(x$tree_row, k = 5))
  return(clust)
}


heatmap_DE_ts <- function(sig, rld, org, tissue, type, title){
  rows <- match(sig$id, row.names(rld))
  mat <- assay(rld)[rows,]
  n=dim(mat)[[1]]
  df <- as.data.frame(colData(rld))
  df = df[-c(2,4,5,6,7)]
  breaksList = seq(-2, 2, by = 0.25)
  if (org=="mouse"){
    #time_colors = list(time=c("0"="grey80", "2"="grey50", "8"="grey30"),condition=c(DSS="grey30", Control="grey70"))
    time_colors = list(time=c("0"="#EFF3FF", "2"="#BDD7E7", "8"="#6BAED6"),condition=c(DSS="grey30", Control="grey70"))
  }
  else{
    #time_colors = list(time=c("0"="grey80", "2"="grey60", "4"="grey40", "5"="grey25"),condition=c(DSS="grey30", Control="grey70"))
    #time_colors = list(condition=c(DSS="grey30", Control="grey70"))
    #time_colors = list(time=c("0"="azure", "2"="cadetblue1", "4"="cornflowerblue", "5"="deepskyblue4"),condition=c(DSS="grey30", Control="grey70"))
    time_colors = list(time=c("0"="#EFF3FF", "2"="#BDD7E7", "4"="#6BAED6", "5"="#2171B5"),condition=c(DSS="grey30", Control="grey70"))
    }
  annot_heatmap = c("DSS1","DSS2","DSS3","Control1","Control2","Control3")
  x <- pheatmap(mat, 
                scale="row",
                cluster_rows=TRUE,
                cluster_cols=FALSE, 
                show_rownames=FALSE,
                #show_colnames=FALSE,
                annotation_col=df,
                #annotation_colors = list(condition=c(DSS="grey30", Control="grey70")),
                annotation_colors = time_colors,
                #annotation_legend = FALSE,
                #annotation_names_col = FALSE,
                #annotation_col=c("DSS","Control"),
                key.title = "Row-wise Z-score",
                #breaks = breaksList,
                border_color = "NA",
                #color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
                #cellwidth=20,
                main=paste(title,",n=",n,sep="") )
  save_pheatmap_pdf(x, paste("figures/heatmap_DE_TS_", org, "_",tissue,"_",type,".pdf", sep=""))
  #,width=9,height=7
  clust <- cbind(mat, cluster = cutree(x$tree_row, k = 6))
  return(clust)
}


plotPCA3D <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE){
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1],
                  PC2 = pca$x[, 2],
                  PC3 = pca$x[, 3],
                  group = group,
                  intgroup.df,
                  name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:3]
    return(d)
  }
  message("Generating plotly plot")
  p <- plotly::plot_ly(data = d,
                       x = ~PC1,
                       y = ~PC2,
                       z = ~PC3,
                       color = group,
                       mode = "markers",
                       type = "scatter3d")
  return(p)
}


return_sig<-function(res1,res2,res3,lfc,pval,padj,biotypes){
  if (!is.null(res1)){
    sig1<-res1[abs(res1$logFC)>lfc & res1$pvalue<pval & res1$padj<padj & res1$biotype %in% biotypes,]
  }
  if (!is.null(res2)){
    sig2<-res2[abs(res2$logFC)>lfc & res2$pvalue<pval & res2$padj<padj & res2$biotype %in% biotypes,]
  }
  if (!is.null(res3)){
    sig3<-res3[abs(res3$logFC)>lfc & res3$pvalue<pval & res3$padj<padj & res3$biotype %in% biotypes,]
  }
  if ( !is.null(res1) & is.null(res2) & is.null(res3) ){
    sig <- sig1
  } 
  else if ( !is.null(res1) & !is.null(res2) & is.null(res3) ){
    sig <-merge(sig1, sig2, by.x='id', by.y='id',all=TRUE)
  }
  else if ( !is.null(res1) & !is.null(res2) & !is.null(res3) ){
    print("here")
    sigtmp <-merge(sig1, sig2, by.x='id', by.y='id', all.x=TRUE,all.y=TRUE)
    sig <-merge(sigtmp, sig3, by.x='id', by.y='id',all.x=TRUE, all.y=TRUE)
  }
  return(sig)
}

return_sig_no_biotype<-function(res1,res2,res3,lfc,pval,padj){
  if (!is.null(res1)){
    sig1<-res1[abs(res1$logFC)>lfc & res1$pvalue<pval & res1$padj<padj ,]
  }
  if (!is.null(res2)){
    sig2<-res2[abs(res2$logFC)>lfc & res2$pvalue<pval & res2$padj<padj ,]
  }
  if (!is.null(res3)){
    sig3<-res3[abs(res3$logFC)>lfc & res3$pvalue<pval & res3$padj<padj,]
  }
  if ( !is.null(res1) & is.null(res2) & is.null(res3) ){
    sig <- sig1
  } 
  else if ( !is.null(res1) & !is.null(res2) & is.null(res3) ){
    sig <-merge(sig1, sig2, by.x='id', by.y='id',all=TRUE)
    sig$logFC=max(sig$logFC.x,sig$logFC.y)
  }
  else if ( !is.null(res1) & !is.null(res2) & !is.null(res3) ){
    print("here")
    sigtmp <-merge(sig1, sig2, by.x='id', by.y='id', all.x=TRUE,all.y=TRUE)
    sig <-merge(sigtmp, sig3, by.x='id', by.y='id',all.x=TRUE, all.y=TRUE)
  }
  return(sig)
}

filter_overlapping_mirs <- function(seqf){
  seqf<-seqf[order(seqf$seq),]
  selected_mirs=c()
  seqf[] <- lapply(seqf, as.character)
  i=2
  while (i<length(seqf$seq)){
    #while (i<10){
    k=1
    a<-grepl(seqf[i-1,2],seqf[i,2])
    if (a==TRUE){
      b=TRUE; k=i;
      while (b==TRUE){
        b<-grepl(seqf[k,2],seqf[k+1,2])
        k=k+1
      }  
      low=i-1; high=k-1;
      if((high-low)>0){
        print(paste(low, high, sep=" "))
        i=k+1
        print(seqf[low:high,])
        sim <- seqf[low:high,] #data frame with similar miRNAs 
        annot_sim <- sim[grep("mir", sim$id), ] # data frame with similar miRNAs annotated in miRBase
        #print(annot_sim)
        if (nrow(annot_sim)==0) { 
          mir <- sim[nrow(sim),]$id
          #print(paste("longest mir, mirdeep found: ",mir,sep=""))
          selected_mirs <- c(selected_mirs, mir)
        }
        else {
          mir <- annot_sim[nrow(annot_sim),]$id
          #print(paste("longest mir from mirbase: ",mir,sep=""))
          selected_mirs <- c(selected_mirs, mir)
        }
      } 
    }
    else{ 
      mir <- seqf[i-1,]$id
      #print(paste("mir not overlapping others, included: ",mir, sep=""))
      selected_mirs <- c(selected_mirs, mir )
      i=i+1
    } 
  }  
  if (i<length(seqf$seq)){
    selected_mirs <- c(selected_mirs, seqf[i,]$id)
  }
  return(selected_mirs)
}  

intersect_results <- function(res1,res2){
  common_mcbs <- as.data.frame(merge(res1, res2, by.x='partial_id', by.y='partial_id'))
  common_mcbs <- common_mcbs[,c(1,2:3,1,9:10,4:8,11:15),] 
  names(common_mcbs)[c(1,4)]<-c("partial_id.x","partial_id.y")
  common_mcbs_seq <- as.data.frame(merge(res1, res2, by.x='seq', by.y='seq'))
  common_mcbs_seq <- common_mcbs_seq[,c(2,1,3,9,1,10,4:8,11:15),] 
  names(common_mcbs_seq)[c(2,5)]<-c("seq.x","seq.y")
  names(common_mcbs_seq)
  names(common_mcbs)
  common_mcbs=rbind(common_mcbs, common_mcbs_seq)
  common_mcbs <- unique(common_mcbs)
  mirs_same_id_seq_1 <- common_mcbs$partial_id.x
  mirs_same_id_seq_2 <- common_mcbs$partial_id.y
  res1=res1[!(res1$partial_id %in% mirs_same_id_seq_1),]
  res2=res2[!(res2$partial_id %in% mirs_same_id_seq_2),]
  common_mcbs_partial_seq <- as.data.frame(merge(res1, res2, by.x='partial_seq', by.y='partial_seq'))
  common_mcbs_partial_seq <- common_mcbs_partial_seq[,c(2,3,1,9,10,1,4:8,11:15),] 
  names(common_mcbs_partial_seq)[c(3,6)]<-c("partial_seq.x","partial_seq.y")
  names(common_mcbs_partial_seq)
  names(common_mcbs)
  common_mcbs=rbind(common_mcbs, common_mcbs_partial_seq)
  return(common_mcbs)
}

get_mean<-function(rld,coldata,org,tissue,type){
  #org="mouse";tissue="colon";rld<-assay(rldc);type="total"  
  m<-rld[,0]
  #m$id=rownames(m)
  if (tissue=="colon"){
    c1<-rownames(coldata[coldata$tissue=="colon" & coldata$condition=="DSS",])
    c2<-rownames(coldata[coldata$tissue=="colon" & coldata$condition=="Control",])
    m<-cbind(m, id=rownames(m), DSS=rowMeans(rld[,c1]), Control=rowMeans(rld[,c2]), Diff=rowMeans(rld[,c1])-rowMeans(rld[,c2]))
  }
  if (tissue=="blood" & org=="mouse"){
    c1<-rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" & coldata$time=="0",])
    c2<-rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" & coldata$time=="2",])
    c3<-rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" & coldata$time=="8",])
    m<-cbind(m, id=rownames(m),DSS_day0=rowMeans(rld[,c1]), DSS_day2=rowMeans(rld[,c2]), DSS_day8=rowMeans(rld[,c3]), Diff=rowMeans(rld[,c3])-rowMeans(rld[,c1]))
  }
  m<-as.matrix(m)
  #m[,c("Diff")]<-as.numeric(levels(m[,c("Diff")]))[m[,c("Diff")]]
  m[,c("Diff")]<-as.numeric(unlist(m[,c("Diff")]))
  m<-as.data.frame(m)
  m$Diff<-as.numeric(levels(m$Diff))[m$Diff]
  write.table(m,file=paste("results/mean_rlog_", org, "_", tissue, "_",type,".tsv",sep=""), sep="\t", row.names=FALSE, quote = F )
  return(m)
}



