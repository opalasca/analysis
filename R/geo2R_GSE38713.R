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

