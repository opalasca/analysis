


data<-get_data("data/mouse_strtie_quant_tpm_reseq.csv", "data/config_mouse.txt")
tpmmr<-data[[1]]; coldatam<-data[[2]]
data<-get_data("data/mouse_strtie_quant_tpm.csv", "data/config_mouse.txt")
tpmm<-data[[1]]; coldatam<-data[[2]]

org="mouse"
tpmmr <- log(tpmmr+1)
tpmmr<-as.data.frame(tpmmr)
tpmm <- log(tpmm+1)
tpmm<-as.data.frame(tpmm)

org="mouse"
tpmm$id=rownames(tpmm)
tpmmr$id=rownames(tpmmr)
tpmm_merged <- as.data.frame(merge(tpmm, tpmmr[c("M1-1","M1-2","M1-3","M4-3","M7-1","M7-2","M7-3","M8-2","M8-3","id")], by.x='id', by.y='id', all.x=TRUE))
tpmm_merged_orth <- as.data.frame(merge(tpmm_merged, mh_orth, by.x='id', by.y='Gene.stable.ID', all.x=TRUE))
write.table(tpmm_merged_orth, file=paste("results/", org, "_log_TPM_plus_reseq.tsv", sep=""), sep=",", row.names=F,quote = FALSE)


#Pig 

data<-get_data("data/pig_strtie_quant_tpm.csv", "data/config_pig.txt")
tpmp<-data[[1]]; coldatap<-data[[2]]

tpm<-tpmp
org="pig"
tpm<-tpm[rowSums( tpm != 0 ) >= 3,] 
tpm <- log(tpm+1)
tpm<-as.data.frame(tpm)
tpm$id=rownames(tpm)
tpm_orth <- as.data.frame(merge(tpm, ph_orth, by.x='id', by.y='Gene.stable.ID', all.x=TRUE))
#write.table(tpm_orth, file=paste("results/", org, "_log_TPM_reseq.tsv", sep=""), sep=",",row.names=T, col.names=NA,quote = FALSE)
write.table(tpm_orth, file=paste("results/", org, "_log_TPM.tsv", sep=""), sep=",", row.names=F,quote = FALSE)

pdf("figures/heatmap_DE_genes_mouse_colon_pc.pdf")
heatmap.2(mat, 
          labCol = annot_heatmap,
          dendrogram="row", 
          Colv=FALSE, trace="none", 
          margin=c(6,9),
          distfun = dist,
          hclustfun = hclust,
          lhei = c(1,4),
          #lwid = c(0.3,3),
          #lmat = rbind(c(0,3),c(2,1),c(0,4)),
          col=brewer.pal(11,"RdBu"),
          #col= colorRampPalette(brewer.pal(8, "Blues"))(25),
          #col=colorRampPalette(c("blue","white","red"))(256),
          scale="row", labRow = FALSE,
          #key.par=list(mar=c(4,4,4,4)),
          key=TRUE, key.title = NA, density.info = "none")
dev.off()