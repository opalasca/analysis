#source("https://bioconductor.org/biocLite.R")
#biocLite("goseq")

library("goseq")
library(GO.db)


go_term_enrichment <- function(res, padjthr, pthr, FCthr, org, tissue){
  res <- res[complete.cases(res),]
  assayed.genes <- res$Human.gene.stable.ID
  de.genes <- res$Human.gene.stable.ID[which( res$padj<pthr & abs(res$logFC)>FCthr) ]
  #de.genes <- as.integer(res$Human.gene.stable.ID[which( res$avg.expr > 20  & abs(res$logFC) > 1.5) ])
  gene.vector <- as.integer(assayed.genes%in%de.genes)
  names(gene.vector)=assayed.genes
  table(gene.vector)
  pwf=nullp(gene.vector,"hg19","ensGene") 
  #head(pwf)
  GO.wall=goseq(pwf,"hg19","ensGene",test.cats=c("GO:MF","GO:BP"))
  GO.wall$padj <- p.adjust(GO.wall$over_represented_pvalue,method="BH")
  GOsig <- subset(GO.wall, padj <= padjthr)
  write.table(GOsig[,c(6,8,4,5)], file=paste("results/go_term_enrichment_", org, "_", tissue,"_", pthr,"_",FCthr, ".tsv", sep=""), sep="\t", row.names = FALSE,quote = FALSE)
  d <- GOsig[1:20,c(6,8,4,5)]
  png(paste("tables/go_top_", org, "_", tissue, ".png", sep=""),height = 30*nrow(d), width =150*ncol(d))
  grid.table(d,rows=NULL)
  dev.off()
  return(GOsig)
}

go_term_enrichment_top_expr <- function(res, padjthr, pthr, FCthr, org, tissue){
  res <- res[complete.cases(res),]
  assayed.genes <- res$Human.gene.stable.ID
  de.genes_excl <- res$Human.gene.stable.ID[which( res$padj<pthr & abs(res$logFC)>FCthr) ]
  #de.genes <- as.integer(res$Human.gene.stable.ID[which( res$avg.expr > 20  & abs(res$logFC) > 1.5) ])
  n=700
  de.genes <- res %>%
    arrange(desc(avg.expr)) %>% 	#Arrange rows by padj values
    pull(Human.gene.stable.ID) %>% 		#Extract character vector of ordered genes
    .[1:n] 	
  de.genes <- as.vector(de.genes)
  de.genes <- de.genes[!de.genes %in% de.genes_excl]
  gene.vector <- as.integer(assayed.genes%in%de.genes)
  names(gene.vector)=assayed.genes
  table(gene.vector)
  pwf=nullp(gene.vector,"hg19","ensGene") 
  #head(pwf)
  GO.wall=goseq(pwf,"hg19","ensGene",test.cats=c("GO:MF","GO:BP"))
  GO.wall$padj <- p.adjust(GO.wall$over_represented_pvalue,method="BH")
  GOsig <- subset(GO.wall, padj <= padjthr)
  #write.table(GOsig[,c(6,8,4,5)], file=paste("results/go_term_enrichment_", org, "_", tissue,"_", pthr,"_",FCthr, ".tsv", sep=""), sep="\t", row.names = FALSE,quote = FALSE)
  d <- GOsig[1:20,c(6,8,4,5)]
  png(paste("tables/go_expr_", org, "_", tissue, ".png", sep=""),height = 30*nrow(d), width =150*ncol(d))
  grid.table(d,rows=NULL)
  dev.off()
  return(GOsig)
}

#enrichment for a given gene list 

go_clust1 <- go_term_enrichment_list(clust,resmb,0.1,"mouse","blood","1")
go_clust2 <- go_term_enrichment_list(clust,resmb,0.1,"mouse","blood","2")

go_clust3 <- go_term_enrichment_list(clust,resmb,1,"mouse","blood","3")
go_clust4 <- go_term_enrichment_list(clust,resmb,0.1,"mouse","blood","4")
go_clust5 <- go_term_enrichment_list(clust,resmb,0.1,"mouse","blood","5")
go_clust6 <- go_term_enrichment_list(clust,resmb,1,"mouse","blood","6")


go_pig_col=go_term_enrichment(respc,0.1,0.2,1, "pig","colon")
go_mouse_col=go_term_enrichment(resmc,0.1,0.1,1,"mouse","colon")
go_pig_blood=go_term_enrichment(respb,0.1,0.1,1,"pig","blood")
go_mouse_blood=go_term_enrichment(resmb,0.1,0.1,1,"mouse","blood")

go_pig_col_expr=go_term_enrichment_top_expr(respc,0.1,0.2,1, "pig","colon")
go_mouse_col_expr=go_term_enrichment_top_expr(resmc,0.1,0.1,1,"mouse","colon")
go_pig_blood_expr=go_term_enrichment_top_expr(respb,0.1,0.1,1,"pig","blood")
go_mouse_blood_expr=go_term_enrichment_top_expr(resmb,0.1,0.1,1,"mouse","blood")

library(gridExtra)
library(grid)

d <- go_pig_blood[1:20,c(6,8,4,5)]
png("tables/go_test.png", height = 30*nrow(d), width =150*ncol(d))
grid.table(d,rows=NULL)
dev.off()

t<-go_pig_blood[,c(6,8,4,5)]
write.table(t, file=paste("results/go_term_enrichment_", org, "_", tissue,"_", pthr,"_",FCthr, ".tsv", sep=""), sep="\t")


enriched.pig=go_pig_col$category[p.adjust(go_pig_col$over_represented_pvalue,method="BH")<.01]
enriched.pig=go_pig_col[p.adjust(go_pig_cnol$over_represented_pvalue,method="BH")<.01]
go_pig_col_adj<-go_pig_col
go_pig_col_adj$padj <- p.adjust(go_pig_col$over_represented_pvalue,method="BH")
enriched.mouse=go_mouse_col$category[p.adjust(go_mouse_col$over_represented_pvalue,method="BH")<.05]
for(go in enriched.pig){
  print(Term(GOTERM[[go]]))
  #print(GOTERM[[go]])
  #cat("--------------------------------------\n") 
}
for(go in enriched.mouse[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n") 
}
enriched.GO=go_mouse_col$category[p.adjust(go_mouse_col$over_represented_pvalue,method="BH")<.05]
enriched.GO

