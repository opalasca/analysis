library(biomaRt)
library(mygene)
library(grid)
library(gridExtra)

source("R/deseq_longRNA.R")


human_id <- (as.vector(resm[resm$id == mouseg,]))[[8]]
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

resm <- resm[complete.cases(resm),]
resp <- resp[complete.cases(resp),]

# For a specific gene, plot the values across samples for mouse and pig 
# and separately human UC set
gene="SLPI"
gene="LGMN"
mouse_gene <- resm$id[resm$Human.gene.name == gene]
pig_gene <- resp$id[resp$Human.gene.name == gene]
if (!identical(pig_gene, character(0)) & !identical(mouse_gene, character(0))){
  plot_two_genes(mouse_gene, gene, ddsm, pig_gene, gene, ddsp)
}
png(paste("figures/", gene, "_human.png", sep=""))
myrow=rownames(annot[annot$Gene.Symbol==gene,])
groups <- gset@phenoData@data[["description"]]
p<-plot(jitter(as.integer(groups)),exprs(gset)[myrow,], ylab="log intensity", xaxt="n",
        xlab="", xlim=c(0.5, nlevels(groups)+0.5), main=gene)
p+axis(side=1, at=seq_len(nlevels(groups)), label=levels(groups))
dev.off()


# Plot top n mouse genes in mouse and pig
resm=resmc
resp=respc
ddsm=ddsmc
ddsp=ddspc
n=40
top_mouse <- resm %>%
  arrange(padj) %>% 	#Arrange rows by padj values
  #arrange(desc(avg.expr)) %>% 	#Arrange rows by padj values
  pull(id) %>% 		#Extract character vector of ordered genes
  .[1:n] 		#Extract the fir??st 20 genes
not_found=0
pdf(file=paste("figures/top_mouse.pdf"),width=8, height=3) 
for (i in 1:n){
  i=2
  mouseg=top_mouse[i]
  pig_gene <- resp$id[resp$Mouse.gene.stable.ID == mouseg]
  gene_name <- (as.vector(resm[resm$id == mouseg,]))[[9]]
  gene_name <- (as.vector(resm[resm$id == mouseg,]))[[9]]
  #human_id <- (as.vector(resm[resm$id == mouseg,]))[[8]]
  #genedesc <- getBM(attributes=c('ensembl_gene_id','entrezgene', 'wikigene_description'), filters = 'ensembl_gene_id', values = c(human_id), mart =ensembl)
  #genedesc$entrezgene
  #if (!identical(pig_gene, character(0))){
  if (!is.na(pig_gene)){
    plot_two_genes(mouseg, gene_name, ddsm, pig_gene, gene_name, ddsp)
  }
  else {
    not_found = not_found + 1
    plot_two_genes(mouseg, gene_name, ddsm, "", gene_name, ddsp)
    #plot_gene(mouseg,gene_name,ddsm,genedesc$entrezgene)
  }
} 
dev.off()
not_found

# Plot read profiles for top n pig genes in mouse and pig 
n=20
top_pig <- resp %>%
  arrange(padj) %>% 	#Arrange rows by padj values
  #arrange(desc(avg.expr)) %>% 	#Arrange rows by padj values
  pull(id) %>% 		#Extract character vector of ordered genes
  .[1:n] 		#Extract the first 20 genes
not_found=0
pdf(file=paste("figures/top_pig.pdf"),width=8, height=3) 
for (i in 1:n){
  pigg=top_pig[i]
  mouse_gene <- resm$id[resm$Pig.gene.stable.ID == pigg]
  gene_name <- (as.vector(resp[resp$id == pigg,]))[[9]]
  if (!is.na(pig_gene)){
    plot_two_genes(mouse_gene, gene_name, ddsm, pigg, gene_name, ddsp)
  }
  else {
    plot_two_genes("", gene_name, ddsm, pigg, gene_name, ddsp)
    not_found = not_found + 1}
}  
dev.off()
not_found

human_id <- (as.vector(resm[resm$id == mouseg,]))[[8]]
human_id <- "54680"
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genedesc <- getBM(attributes=c('ensembl_gene_id','entrezgene', 'wikigene_description'), filters = 'entrezgene', values = c(human_id), mart =ensembl)
genedesc$ensembl_gene_id

library("ggpubr")
#test
#mouseg="ENSMUSG00000022385"
#gene_symbol="REG1B"
#groups <- gset@phenoData@data[["description"]]
#myrow=rownames(annot[annot$Gene.Symbol==gene_symbol,])
#array_expr <- exprs(gset)[myrow,]
#pig_gene <- resp$id[resp$Mouse.gene.stable.ID == mouseg]
#pig_gene="ENSSSCG00000005243"
#pig_gene="ENSSSCG00000005479"
#plot_three_genes("ENSSSCG00000005479","IL10",ddsp,"ENSSSCG00000005479","IL10",ddsp,array_expr, groups)
#plot_array_int("bla", array_expr, groups)

# Plot individual genes #
#########################
gene="SLPI"
#gene="LGMN"
mouse_gene <- resm$id[resm$Human.gene.name == gene]
pig_gene <- resp$id[resp$Human.gene.name == gene]
if (!identical(pig_gene, character(0)) & !identical(mouse_gene, character(0))){
  plot_two_genes(mouse_gene, gene, ddsm, pig_gene, gene, ddsp)
}
png(paste("figures/", gene, "_human.png", sep=""))
myrow=rownames(annot[annot$Gene.Symbol==gene,])
groups <- gset@phenoData@data[["description"]]
p<-plot(jitter(as.integer(groups)),exprs(gset)[myrow,], ylab="log intensity", xaxt="n",
        xlab="", xlim=c(0.5, nlevels(groups)+0.5), main=gene)
p+axis(side=1, at=seq_len(nlevels(groups)), label=levels(groups))
dev.off()

# Plot top genes #
##################
resh <- AI_vs_C
top_human <- resh %>%
  #arrange(desc(B)) %>% 	#Arrange rows by padj values
  arrange(logFC) %>% 	#Arrange rows by padj values
  pull(Gene.Symbol) %>% 		#Extract character vector of ordered genes
  .[1:30] 		#Extract the fir??st 20 genes
#pdf(file=paste("figures/top_human.pdf"),width=8, height=3) 
for (i in 1:20){
  gene=top_human[i]
  if (grepl("/",gene)==FALSE){
    print(gene)
    mouse_gene <- resm$id[resm$Human.gene.name == gene]
    pig_gene <- resp$id[resp$Human.gene.name == gene]
    if (!identical(pig_gene, character(0)) & !identical(mouse_gene, character(0))){
      plot_two_genes(mouse_gene, gene, ddsm, pig_gene, gene, ddsp)
    }
    png(paste("figures/", gene, "_human.png", sep=""))
    myrow=rownames(resh[resh$Gene.Symbol==gene,])
    groups <- gset@phenoData@data[["description"]]
    p<-plot(jitter(as.integer(groups)),exprs(gset)[myrow,], ylab="log intensity", xaxt="n",
            xlab="", xlim=c(0.5, nlevels(groups)+0.5), main=gene)
    p+axis(side=1, at=seq_len(nlevels(groups)), label=levels(groups))
    dev.off()
  }
}
dev.off()

resh <- AI_vs_C
top_human <- resh %>%
  arrange(abs(logFC)) %>% 	#Arrange rows by padj values
  pull(Gene.Symbol)
for (i in 1:len(top_human)){
  gene=top_human[i]
  if (grepl("/",gene)==FALSE){
    print(gene)
    mouse_gene <- resm$id[resm$Human.gene.name == gene]
    pig_gene <- resp$id[resp$Human.gene.name == gene]
    if (!identical(pig_gene, character(0)) & !identical(mouse_gene, character(0))){
      plot_two_genes(mouse_gene, gene, ddsm, pig_gene, gene, ddsp)
    }
    png(paste("figures/", gene, "_human.png", sep=""))
    myrow=rownames(resh[resh$Gene.Symbol==gene,])
    groups <- gset@phenoData@data[["description"]]
    p<-plot(jitter(as.integer(groups)),exprs(gset)[myrow,], ylab="log intensity", xaxt="n",
            xlab="", xlim=c(0.5, nlevels(groups)+0.5), main=gene)
    p+axis(side=1, at=seq_len(nlevels(groups)), label=levels(groups))
    dev.off()
  }
}

