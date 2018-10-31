# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Sun Apr 15 08:55:48 EDT 2018
# http://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day1/bioc-intro.pdf

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)
library("dplyr")

# load series and platform data from GEO
acc="GSE38713"
gset <- getGEO("GSE38713", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "0000000000000111111112323232323232322222222"
#            G0           G1  G2 G3
sml <- c()
labels <- c("control", "remission", "active involved", "active_noninvolved")
conditions <- c("C","R","AI","AN")
for (i in 1:nchar(gsms)) { sml[i] <- conditions[as.integer(substr(gsms,i,i))+1] }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

annot <- fData(gset)
annot <- annot[annot$SPOT_ID!="--Control",]

#remove control probes and other probes not of interest (only keep human related probesets)
fData(gset)<-fData(gset)[fData(gset)$SPOT_ID!="--Control",]
rows_to_keep<-rownames(fData(gset))
gset<-gset[rows_to_keep,]

fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)

o <- order(fit$Amean, decreasing=TRUE)
dup <- duplicated(fit$genes$Gene.Symbol[o])
fit.unique <- fit[o,][!dup,]

cont.matrix <- makeContrasts(AI-C, AN-C, AI-AN, AI-R, AN-R, levels=design)
cont_names<-c("AI_vs_C", "AN_vs_C", "AI_vs_AN", "AI_vs_R", "AN_vs_R")
fit2 <- contrasts.fit(fit.unique, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

#pthr = 0.1; logFCthr = 1.5
pthr = 1; logFCthr = 0

for (i in 1:length(cont_names)) {
  tname=cont_names[i]
  tname_merged=paste(cont_names[i],"merged",sep="_")
  t=data.frame()
  t=topTable(fit2, coef=i, adjust="fdr", sort.by="B", number=dim(fit2$coefficients)[1])
  #t<-t[t$adj.P.Val<thr,]
  t<-t[order(t$adj.P.Val),]
  names(t)[names(t) == "adj.P.Val"] = "padj"
  names(t)[names(t) == "AveExpr"] = "avg.expr"
  t<-t[t$padj < pthr,]
  t<-t[abs(t$logFC) > logFCthr,]
  t<- t[c(11,17:22,1,2,10,12)]
  assign(tname,t)
  write.table(t, file=paste("results/",acc,"_",tname, ".tsv",sep=""), row.names=F, sep="\t")
}

#library(biomaRt)
#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#genedesc <- getBM(attributes=c('ensembl_gene_id', 'entrezgene', 'wikigene_description'), filters = 'entrezgene', values = c(human_id), mart =ensembl)
#genedesc$ensembl_gene_id

#Transform avg intensities per group to percentiles 
avg_intensities <- as.data.frame(fit.unique[["coefficients"]])
avg_intensities$mean <- fit.unique$Amean 
avg_intensities <- merge(avg_intensities,annot["Gene.Symbol"],by="row.names",all.x=TRUE)

data<-get_data("data/pig_strtie_quant_tpm.csv", "data/config_pig.txt")
tpmp<-data[[1]]; coldatap<-data[[2]]
data<-get_data("data/mouse_strtie_quant_tpm.csv", "data/config_mouse.txt")
tpmm<-data[[1]]; coldatam<-data[[2]]

tpmm <- transform_tpm(tpmm, coldatam, "mouse") 
tpmp <- transform_tpm(tpmp, coldatap, "pig") 

getDE <- function(tpm, samp1, samp2){
  #tpm<-tpmp;samp1="P20-1"; samp2=
  res <- tpm[,c("id", samp1, samp2, "Human.gene.stable.ID", "Human.gene.name")]
  res$logFC <- as.numeric(tpm$samp2) / as.numeric(tpm$samp1)
  return(res)
} 

resP20 <- getDE(tpmp,"P20-1","P20-4")
#avg_intensities <- as.data.frame(fit.unique[["coefficients"]])
#avg_intensities$mean <- fit.unique$Amean 

IPA <- read.csv('data/IBD_IPA.txt', sep="\t")
IBD_Reza <- read.csv('data/IBD_genes_Reza.txt', sep="\t", header=FALSE)
IBD_GWAS <- read.csv('data/IBD_associated_loci_Reza.txt', sep="\t", header=FALSE)
common_hm <- data.frame()
common_hp <- data.frame()
common_mp <- data.frame()
#res1<-as.data.frame(merge(avg_intensities, AI_vs_C, by.x='Gene.Symbol', by.y='Gene.Symbol'))
res1<-avg_intensities
res1<-as.data.frame(merge(avg_intensities, IPA, by.x='Gene.Symbol', by.y='Molecule.Name'))
res2<-tpmm
res3<-tpmp
common_hm <- as.data.frame(merge(res1, res2, by.x='Gene.Symbol', by.y='Human.gene.name'))
common_hp <- as.data.frame(merge(res1, res3, by.x='Gene.Symbol', by.y='Human.gene.name'))
common_mp <- as.data.frame(merge(res2, res3, by.x='Human.gene.name', by.y='Human.gene.name'))
common_mp_filtered <- as.data.frame(merge(common_mp, resm_filtered, by.x='Human.gene.name', by.y='Human.gene.name'))

get_correlations(common_hm, "human", "mouse") 
get_correlations(common_hp, "human", "pig") 

common_mp = common_mp_filtered
cor(common_mp$qCD.x, common_mp$qCD.y)
cor(common_mp$qCD.x, common_mp$qCC.y)
cor(common_mp$qCC.x, common_mp$qCC.y)
cor(common_mp$qCC.x, common_mp$qCD.y)

require("lsa")
cosine(as.matrix(common_hm[,c("AI","C","CD","CC")]))
cor(as.matrix(common_hm[,c("AI","C","CD","CC")]))
cosine(as.matrix(common_hp[,c("AI","C","CD","CC")]))
cor(as.matrix(common_hp[,c("AI","C","CD","CC")]))
cosine(as.matrix(common_mp[,c("CD.x","CC.x","CD.y","CC.y")]))
cor(as.matrix(common_mp[,c("CD.x","CC.x","CD.y","CC.y")]))

cosine(as.matrix(common_mp[,c("BD.x","BC.x","BD.y","BC.y")]))
cor(as.matrix(common_mp[,c("BD.x","BC.x","BD.y","BC.y")]))

drops <- c("Human.gene.stable.ID.x","Human.gene.stable.ID.y","Human.gene.name","id.x","id.y")
cos_mouse_pig <- cosine(as.matrix(common_mp[, !(names(common_mp) %in% drops)]))
write.table(cos_mouse_pig, file="results/mouse_pig_cosine.tsv", row.names=T,col.names=NA, sep="\t")



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

png(paste("figures/hist_AI.png", sep=""))
qplot(AI, data=avg_intensities, geom='histogram')
dev.off()
png(paste("figures/hist_C.png", sep=""))
qplot(C, data=avg_intensities, geom='histogram')
dev.off()

png(paste("figures/hist_avg_col_mouse.png", sep=""))
qplot(CD, data=tpmm, geom='histogram')
dev.off()
png(paste("figures/hist_avg_col_pig.png", sep=""))
qplot(CD, data=tpmp, geom='histogram')
dev.off()

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
