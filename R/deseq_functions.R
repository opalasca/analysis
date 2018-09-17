#install.packages("Vennerable", repos="http://R-Forge.R-project.org") 
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
  all(rownames(coldata) %in% colnames(cts))
  cts <- cts[, rownames(coldata)]
  all(rownames(coldata) == colnames(cts))
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

get_deseq <- function(cts, coldata, sample_list, rowsums){
  dds <- DESeqDataSetFromMatrix(countData = cts[,sample_list],
                              colData = coldata[sample_list,],
                              design = ~ condition)
                              #design = ~ time)
  print(nrow(dds))
  dds<-dds[rowSums( counts(dds) != 0 ) >= rowsums,] 
  print(nrow(dds))
  dds$condition <- relevel(dds$condition, ref = "C")
  dds <- DESeq(dds)
  print(nrow(dds))
  return(dds)
}

get_deseq_blood <- function(cts, coldata, from, to, rowsums){
  dds <- DESeqDataSetFromMatrix(countData = cts[,c(from:to)],
                                colData = coldata[from:to,],
                                #design = ~ condition + time + condition:time)
                                design = ~ time + condition + time:condition)
  print(nrow(dds))
  dds<-dds[rowSums( counts(dds) != 0 ) >= rowsums,] 
  print(nrow(dds))
  dds$condition <- relevel(dds$condition, ref = "C")
  #dds <- DESeq(dds, full= ~ condition + time + condition:time, reduced = ~time + condition, test="LRT")
  dds <- DESeq(dds, test="LRT", reduced= ~ time + condition)
  return(dds)
}

# draw basic data exploration plots - heatmap, sample distance, PCA. 
# vsd can be replaced with other normalization results - e.g. rld, ntd   
basic_plots <-function(dds, vsd, org, tissue, type){ 
  #Heatmap
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:20]
  df <- as.data.frame(colData(dds)[,c("condition","time")])
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
 
  plotPCA(vsd, intgroup=c("condition")) +
    #geom_point(aes(size=25))+
    aes(label = name)+
    geom_text_repel()+
    theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"),
          plot.title = element_text(hjust = 0.5,size=20,face="bold"),legend.text = element_text(size=12, face="bold"))
  ggsave(paste("figures/PCA_", org, "_",tissue,"_",type, ".png", sep=""), width=5, height=3.5, dpi = 1000)
  
}

get_results <- function(dds, org, tissue, type, thr){
  res_uc_co <- results(dds)#, contrast=c("condition","D","C"))
  resOrdered_uc_co <- as.data.frame(res_uc_co[order(res_uc_co$pvalue),])
  names(resOrdered_uc_co)[names(resOrdered_uc_co) == "log2FoldChange"] = "logFC"
  names(resOrdered_uc_co)[names(resOrdered_uc_co) == "baseMean"] = "avg.expr"
  #resOrdered_uc_co <- subset(resOrdered_uc_co, padj<thr)
  #resOrdered_uc_co <- subset(resOrdered_uc_co, abs(logFC)>2) 
  resOrdered_uc_co$common_id <-tolower(gsub("_.*","",rownames(resOrdered_uc_co)))
  write.table(resOrdered_uc_co[c(1,2,6,7)], file=paste("results/", org, "_", tissue, "_",type, "_", 1, "_thr", ".tsv", sep=""), sep="\t")
  resOrdered_uc_co <- subset(resOrdered_uc_co, padj<thr)
  return(resOrdered_uc_co)
  
}

pca_plot <- function(vsd, org, tissue, type){
  #org="mouse"; tissue="colon"; type="small"
  plotPCA(vsd, intgroup=c("condition")) + 
    aes(label = name,size=9)+
    geom_text_repel()+
    theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),plot.title = element_text(hjust = 0.5,size=20,face="bold"))
  ggsave(paste("figures/PCA_", org, "_",tissue,"_",type, ".png", sep=""), width=9, height=7, dpi = 100)

}

plot_gene <- function(gene,dds){
  d <- plotCounts(dds, gene=gene, intgroup="condition", returnData=TRUE)
  # Plotting the gene normalized counts, using the samplenames (rownames of d as labels)
  ggplot(d, aes(x = condition, y = count, color = condition)) + 
    geom_point(position=position_jitter(w = 0.1,h = 0)) +
    geom_text_repel(aes(label = rownames(d))) + 
    theme_bw() +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5))
} 

get_orth_table <- function(orth_df, test_col, keep_cols){
  orth_o2o <- orth_df[orth_df[,test_col]=="ortholog_one2one",][,keep_cols]
  orth_o2o <- orth_o2o[!duplicated(orth_o2o), ]
  return(orth_o2o)
}

