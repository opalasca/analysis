#source("https://bioconductor.org/biocLite.R")
#biocLite("tximportData")
#source("https://bioconductor.org/biocLite.R")
#biocLite("ensembldb")
#biocLite("AnnotationHub")

setwd("~/Desktop/IBD/BGI_analysis")


library(AnnotationHub)
ah <- AnnotationHub()
query(ah, "EnsDb.SScrofa")
edb <- ah[["AH57796"]]
edb <- ah[["AH56713"]]
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
tx2gene <- transcripts(edb,return.type="DataFrame")[c("tx_id","gene_id")]

tx2gene <- read.csv("data/tx2gene.SScrofa11.1_new",sep="\t",header=FALSE)
colnames(tx2gene) <- c("tx_id","gene_id")
tx2gene

library(tximportData)
#dir <- system.file("extdata", package = "tximportData")
#list.files(dir)

config_file = "data/config_pig.txt"
dir="data/salmon_pig"
samples <- read.table("data/config_pig.txt", header = TRUE)
coldata <- read.csv(config_file, sep="\t", row.names=1)
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
coldata$time <- as.factor(coldata$time)
sample_list <- colon_samples

files <- file.path(dir, sample_list, "quant.sf.gz")
names(files) <- colon_samples
all(file.exists(files))

library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = coldata[sample_list,],
                                   design = ~ condition)
print(nrow(ddsTxi))
ddsTxi<-ddsTxi[rowSums( counts(ddsTxi) != 0 ) >= 1,] 
print(nrow(ddsTxi))
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "C")
dds <- DESeq(ddsTxi)
res <- results(dds)
res <- results(dds, independentFiltering = FALSE )
res$id <- rownames(res)
res <- as.data.frame(res)
res$id <- gsub("\\..*","",res$id)

resp2 <- as.data.frame(get_results(dds, "pig", "colon", "total", 1, ph_orth, pm_orth))

ddsp<-dds
n=10
top_pig <- resp2[resp2$avg.expr>50,] %>%
  arrange(desc(logFC)) %>% 	#Arrange rows by padj values
  #arrange(desc(avg.expr)) %>% 	#Arrange rows by padj values
  pull(id) %>% 		#Extract character vector of ordered genes
  .[1:n] 		#Extract the first 20 genes
pdf(file=paste("figures/top_pig_salmon.pdf"),width=8, height=3) 

for (i in 1:n){
  pigg=top_pig[i]
  gene_name="none"
  plot_two_genes("", gene_name, ddsm, pigg, gene_name, ddsp)
}  
dev.off()



res$gene <- rownames(res)
res <- as.data.frame(res[order(res$padj),])
plotCounts(dds, gene="ENSSSCG00000035520.1", intgroup="condition")
ggsave("figures/testsalmon.png", width=5, height=4, dpi = 1000)
plot_two_genes("", "test", dds, "ENSSSCG00000012600.3", "test", dds)

print(nrow(dds))
ddsp <- dds
vsd <- vst(ddsp, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
meanSdPlot(assay(vsd))
basic_plots(ddsp, vsd, "pig", "colon","total")
resp <- as.data.frame(get_results(ddsp, "pig", "colon", "total", 1, ph_orth, pm_orth))
