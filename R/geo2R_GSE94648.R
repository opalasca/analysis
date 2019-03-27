# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Fri Jan 4 09:37:42 EST 2019

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)
library(biomaRt)

# load series and platform data from GEO
# Human long RNA array blood
acc="GSE94648"
gset <- getGEO("GSE94648", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL19109", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

 #group names for all samples
gsms <- paste0("11111111111111111222222220000000000000000000000333",
               "33333333333333333333333333333333333333444444444")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
labels <- c("Control","aUC","iUC","aCD","iCD")
conditions <- c("C","AUC","IUC","ACD","ICD")
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
# use file downloaded from ncbi with entrez gene - symbol correspondence 
# cat gene_result.txt  | cut -f3,6 > hsa_gene_id_name.txt

fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)

o <- order(fit$Amean, decreasing=TRUE)
dup <- duplicated(fit$genes$ENTREZ_GENE_ID[o])
fit.unique <- fit[o,][!dup,]

cont.matrix <- makeContrasts(AUC-C, IUC-C, ACD-C, ICD-C, levels=design)
cont_names<-c("AUC_vs_C", "IUC_vs_C", "ACD_vs_C", "ICD_vs_C")
fit2 <- contrasts.fit(fit.unique, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

pthr = 0.1; logFCthr = 1
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
  #t<- t[c(11,17:22,1,2,10,12)]
  assign(tname,t)
  write.table(t, file=paste("results/",acc,"_",tname, ".tsv",sep=""), row.names=F, sep="\t")
}

library(biomaRt)
human_id="2207"
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genedesc <- getBM(attributes=c('ensembl_gene_id', 'entrezgene', 'symbol','wikigene_description'), filters = 'entrezgene', values = c(human_id), mart =ensembl)
#genedesc$ensembl_gene_id


