source("https://bioconductor.org/biocLite.R")
biocLite("goseq")

library("goseq")

#dds<-ddsp
#res <- results(dds, independentFiltering = FALSE )
res <- resp
res <- res[complete.cases(res),]
 
assayed.genes <- res$Human.gene.stable.ID
de.genes <- as.integer(res$Human.gene.stable.ID[which( res$avg.expr > 20  & abs(res$logFC) > 1.5) ])
assayed.genes <- res$Human.gene.stable.ID

#assayed.genes <- res$Human.gene.stable.ID
genes <- as.integer(which( res$avg.expr > 20  & abs(res$logFC) > 1.5) )
names(genes) <- res$Human.gene.stable.ID

pwf=nullp(genes,"hg19","ensGene") 
head(pwf)