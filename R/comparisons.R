
source("../integrated_miRNA_analysis/R/geo2R_GSE48957.R")
source("../integrated_miRNA_analysis/R/geo2R_GSE32273.R")
source("R/geo2R_GSE94648.R")
source("R/geo2R_GSE38713.R")

setwd("~/Desktop/IBD/BGI_analysis")

curated_mirs_blood <- read.csv("data/ranked_Blood_UC_vs_CO.tsv", header=FALSE, sep="\t")
curated_mirs_blood <-curated_mirs_blood[-c(5,3)]
names(curated_mirs_blood) <- c("occ","partial_id","direction")
curated_mirs_blood$partial_id=tolower(curated_mirs_blood$partial_id)

#Mouse colon vs blood total
res2<-resmc
res3<-resmb8
common_mcbt <- data.frame()
common_mcbt <- as.data.frame(merge(res2, res3, by.x='id', by.y='id'))
common_mcbt <- common_mcbt[(common_mcbt$biotype.x=="protein_coding"),]
common_mcbt <- common_mcbt[(common_mcbt$biotype.x %in% lncRNAs),]
pthr=0.1; thr=1
mp<-get_stats(common_mcbt,  "colon", "blood", pthr, pthr, thr, thr, 0)
mp[[1]]
cons_mcbt<-mp[[2]]
dim(cons_mcbt[cons_mcbt$logFC.x<1,])
cons_lnc <- cons_mcbt[cons_mcbt$biotype.x %in% lncRNAs,]
cons_pc <- cons_mcbt[cons_mcbt$biotype.x=="protein_coding",]


#Pig colon vs blood total
res2<-respc
res3<-respb4
common_mcbt <- data.frame()
common_mcbt <- as.data.frame(merge(res2, res3, by.x='id', by.y='id'))
common_mcbt <- common_mcbt[(common_mcbt$biotype.x=="protein_coding"),]
#common_mcbt <- common_mcbt[(common_mcbt$biotype.x %in% lncRNAs),]
pthr=0.1; thr=1
mp<-get_stats(common_mcbt,  "colon", "blood", pthr, pthr, thr, thr, 0)
mp[[1]]
cons_mcbt<-mp[[2]]

library(knitr) 
mph<-kable(mp[[1]], format = "html")
library(htmlwidgets)
saveWidget(mph, file="results/m.html")

#Mouse colon vs blood small
res2<-resmcs
res3<-resmbs
common_mcbs <- data.frame()
common_mcbs <- as.data.frame(merge(res2, res3, by.x='id', by.y='id'))
pthr=0.05; thr=0
mp<-get_stats_pval(common_mcbs, "colon", "blood", pthr, pthr, thr, thr, 0)
mp[[1]]
cons_mcbs<-mp[[2]]

#Mouse blood day2 vs blood day 8 small
res2<-resmbs2
res3<-resmbs8
common_mcbs <- data.frame()
common_mcbs <- as.data.frame(merge(res2, res3, by.x='id', by.y='id'))
pthr=0.1; thr=0
mp<-get_stats(common_mcbs, "blood2", "blood8", pthr, pthr, thr, thr, 0)
kable(mp[[1]])
cons_mcbs<-mp[[2]]

#Mouse colon vs human UC colon GSE89667
res2<-resmcs
res3<-resOrdered_UC_DD
common_mcbt <- data.frame()
common_mcbt <- as.data.frame(merge(res2, res3, by.x='partial_id', by.y='partial_id'))
pthr=0.1; thr=1
mp<-get_stats(common_mcbt,  "mouse", "human", pthr, pthr, thr, thr, 0)
kable(mp[[1]])
cons_mcbt<-mp[[2]]

#Mouse colon vs human UC colon GSE48957
res2<-resmcs
res3<-UCa_vs_C_merged
common_mcbt <- data.frame()
common_mcbt <- as.data.frame(merge(res2, res3, by.x='seq', by.y='Sequence'))
pthr=0.1; thr=0
mp<-get_stats(common_mcbt,  "mouse", "human", pthr, pthr, thr, thr, 0)
mp[[1]]
cons_mcbt<-mp[[2]]

#Mouse colon vs mouse blood small
res2<-resmcs
res3<-resmbs
common_mcbt <- data.frame()
common_mcbt <- as.data.frame(merge(res2, res3, by.x='partial_id', by.y='partial_id'))
pthr=0.05; thr=0
mp<-get_stats_pval(common_mcbt,  "colon", "blood", pthr, pthr, thr, thr, 0)
kable(mp[[1]])
cons_mcbt<-mp[[2]]

# LogFC comparisons between mouse and pig 

#Mouse vs pig total colon
res2<-resmc[complete.cases(resmc),]
res3<-respc[complete.cases(respc),]
common_mp <- data.frame()
common_mp <- as.data.frame(merge(res2, res3, by.x='Human.gene.name', by.y='Human.gene.name'))
pthr=0.1; thr=0.5
mp<-get_stats(common_mp, "mouse", "pig", pthr, 0.1, thr, thr, 0)
mp<-get_stats(common_mp, "mouse", "pig", 0.1, 0.1, thr, thr, 0)
mp[[1]]
cons_mpc<-mp[[2]]




# Human vs mouse vs pig 
res1<-AI_vs_C
res2<-resmc[complete.cases(resmc),]
res3<-respc[complete.cases(respc),]
#res1<-as.data.frame(merge(AI_vs_C, IPA, by.x='Gene.Symbol', by.y='Molecule.Name'))
#res2<-as.data.frame(merge(res2, IPA, by.x='Human.gene.name', by.y='Molecule.Name'))
#res3<-as.data.frame(merge(res3, IPA, by.x='Human.gene.name', by.y='Molecule.Name'))
common_hm <- data.frame()
common_hp <- data.frame()
common_mp <- data.frame()
common_hm <- as.data.frame(merge(res1, res2, by.x='Gene.Symbol', by.y='Human.gene.name'))
common_hp <- as.data.frame(merge(res1, res3, by.x='Gene.Symbol', by.y='Human.gene.name'))
common_mp <- as.data.frame(merge(res2, res3, by.x='Human.gene.name', by.y='Human.gene.name'))
pthr=0.1; thr=1
hm<-get_stats(common_hm, "human", "mouse", pthr, pthr, thr, thr,0)
hp<-get_stats(common_hp, "human", "pig", pthr, 0.1, thr, thr, 0)
mp<-get_stats(common_mp, "mouse", "pig", pthr, 0.1, thr, 1, 0)

cons_hm<-hm[[2]]
cons_hp<-hp[[2]]
cons_mp<-mp[[2]]
hm[[1]]
hp[[1]]
mp[[1]]



#Mouse vs pig total blood
res2<-resmb[complete.cases(resmb),]
res3<-respb[complete.cases(respb),]
#res3<-respb5[complete.cases(respb5),]
common_mp <- data.frame()
common_mp <- as.data.frame(merge(res2, res3, by.x='Human.gene.name', by.y='Human.gene.name'))
pthr=0.1; thr=1
mp<-get_stats(common_mp, "mouse", "pig", pthr, 0.1, thr, thr, 0)
kable(mp[[1]])
cons_mpb<-mp[[2]]



pthr=0.1; thr=0
#Mouse blood small vs human UC blood GSE32273
res1<-resmbs
#res2<-UCpl_vs_Cpl_merged
#res2<-UCP_vs_CP_merged
res2<-UCm_vs_Cm_merged
names(res2)[10]="seq"
res2$seq=toupper(res2$seq)
res2$seq<-gsub("T", "U", res2$seq)
res1$partial_seq=substr(res1$seq,1,20)
res2$partial_seq=substr(res2$seq,1,20)
res1<-res1[,c(8,10,11,1:3,6:7)]
res2<-res2[,c(13,10,14,1,6,7,9,8)]
common_mcbs <- intersect_results(res1,res2)
mp<-get_stats_pval(common_mcbs,  "mouse", "human", pthr,pthr,thr,thr, 0)
kable(mp[[1]])
cons_mhbs<-mp[[2]]
kable(mp[[2]][c(1,4,9,14)])
kable(mp[[3]][c(1,4,9,14)])

hc<-merge(res1, curated_mirs_blood, by="partial_id")
res=res1
i=1
n=0
for (i in 1:nrow(curated_mirs_blood)){
  id=curated_mirs_blood[i,]$partial_id
  #print(res[i,]$id)
  #print(seq)
  sim <- res[grep(id,res$partial_id), ] 
  #sim <- sim[c(1,5:7)]
  sim <- sim[c(8,1:3,6)]
  if (nrow(sim[sim$pvalue<0.1,])>0) {
    print(curated_mirs_blood[i,])
    print(sim[sim$pvalue<0.1,])
    n=n+1
    #names(res2[i,])<- paste( names(res2[i,]), ".x", sep="")
  }
}  
print(n)


# Pig blood small vs human UC blood GSE32273
res1=respbs4corr
#res2<-UCpl_vs_Cpl_merged
res2<-UCP_vs_CP_merged
res2<-UCm_vs_Cm_merged
names(res2)[10]="seq"
res2$seq=toupper(res2$seq)
res2$seq<-gsub("T", "U", res2$seq)
res1$partial_seq=substr(res1$seq,1,20)
res2$partial_seq=substr(res2$seq,1,20)
res1<-res1[,c(8,10,11,1:3,6:7)]
res2<-res2[,c(13,10,14,1,6,7,9,8)]
common_mcbs <- intersect_results(res1,res2)
mp<-get_stats_pval(common_mcbs,  "pig", "human", pthr, pthr, thr, thr, 0)
kable(mp[[1]])
cons_phbs<-mp[[2]]
kable(mp[[2]][c(1,4,9,14)])
kable(mp[[3]][c(1,4,9,14)])

#Mouse vs pig blood small
res1<-resmbs
res2<-respbs4corr
res1$partial_seq=substr(res1$seq,1,20)
res2$partial_seq=substr(res2$seq,1,20)
res1<-res1[,c(8,10,11,1:3,6:7)]
res2<-res2[,c(8,10,11,1:3,6:7)]
common_mcbs <- intersect_results(res1,res2)
#pthr=0.1; thr=0
mp<-get_stats_pval(common_mcbs, "mouse", "pig", pthr, pthr, thr, thr, 0)
kable(mp[[1]])
cons_mpbs<-mp[[2]]
kable(mp[[2]][c(1,4,9,14)])
kable(mp[[3]][c(1,4,9,14)])

overlap1<-merge(cons_mpbs,cons_mhbs,by="partial_id.x")
overlap2<-merge(cons_mpbs,cons_phbs, by.x="partial_id.y", by.y="partial_id.x")
overlap3<-merge(cons_mhbs,cons_phbs, by.x="partial_id.y", by.y="partial_id.y")


pthr=0.1; thr=0
#Mouse vs pig colon small
res1<-resmcs
res2<-respcscorr
res1$partial_seq=substr(res1$seq,1,20)
res2$partial_seq=substr(res2$seq,1,20)
res1<-res1[,c(8,10,11,1:3,6:7)]
res2<-res2[,c(8,10,11,1:3,6:7)]
common_mccs <- intersect_results(res1,res2)
#pthr=0.1; thr=0
mp<-get_stats_pval(common_mccs, "mouse", "pig", pthr, pthr, thr, thr, 0)
kable(mp[[1]])
cons_mpcs<-mp[[2]]
kable(mp[[2]][c(1,4,9,14)])
kable(mp[[3]][c(1,4,9,14)])  



#selected <- subset(common_hp, abs(common_hp$logFC.x)>2 & abs(common_hp$logFC.y)>2) #& common$logFC.x*common$logFC.y > 0
cor(common_hm$logFC.x, common_hm$logFC.y, method="spearman")
cor(common_hp$logFC.x, common_hp$logFC.y, method="spearman")
cor(common_mp$logFC.x, common_mp$logFC.y, method="spearman")

draw_logFC_corr("test_human", "mouse_IPA", common_hm, 0.05, 0.05)

draw_logFC_corr("human", "mouse_IPA", common_hm, 0.05, 0.05)
draw_logFC_corr("human", "pig_IPA", common_hp, 0.05, 0.05)
draw_logFC_corr("mouse", "pig_IPA", common_mp, 0.05, 0.05)

draw_logFC_corr("human", "mouse", common_hm, 0.05, 0.05)
draw_logFC_corr("human", "pig", common_hp, 0.05, 0.05)
draw_logFC_corr("mouse", "pig", common_mp, 0.05, 0.05)

#use full intersection file for plotting the disjoint sets as well
#common_full <- as.data.frame(merge(res1, res2, by.x=common1, by.y=common2, all.x=TRUE,all.y=TRUE))
#rownames(common)=common$common_id
write.table(common, file=paste("results/common_",out_name,".tsv",sep=""), sep="\t")
draw_logFC_corr(acc1, acc2, cp1, cp2, common, 0.1, 0.1)
draw_volcano_highlight_uncommon(acc1, acc2, cp1, cp2, common_full, res1, res2, 0.1, 0.1)

i=1
n=0
for (i in 1:nrow(res2)){
  seq=substr(res2[i,]$seq,1,20)
  #print(res[i,]$id)
  #print(seq)
  sim <- res3[grep(seq, res3$seq), ] 
  if (nrow(sim)>0) {
    print(res[i,]$id)
    print(sim)
    n=n+1
    names(res2[i,])<- paste( names(res2[i,]), ".x", sep="")
  }
}  
print(n)


