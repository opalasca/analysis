source("../integrated_miRNA_analysis/R/geo2R_GSE48957.R")
source("../integrated_miRNA_analysis/R/geo2R_GSE32273.R")



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
mp[[1]]
cons_mcbs<-mp[[2]]

#Mouse colon vs human UC colon GSE89667
res2<-resmcs
res3<-resOrdered_UC_DD
common_mcbt <- data.frame()
common_mcbt <- as.data.frame(merge(res2, res3, by.x='partial_id', by.y='partial_id'))
pthr=0.1; thr=1
mp<-get_stats(common_mcbt,  "mouse", "human", pthr, pthr, thr, thr, 0)
mp[[1]]
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

#Mouse blood vs human UC blood GSE32273
res2<-resmbs
#res3<-UCpl_vs_Cpl_merged
res3<-UCm_vs_Cm_merged
#res3<-UCP_vs_CP_merged
common_mcbt <- data.frame()
common_mcbt <- as.data.frame(merge(res2, res3, by.x='partial_id', by.y='partial_id'))
pthr=0.3; thr=0
mp<-get_stats(common_mcbt,  "mouse", "human", pthr, pthr, thr, thr, 0)
mp[[1]]
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

#Mouse vs pig total blood

res2<-resmb8[complete.cases(resmb8),]
res3<-respb4[complete.cases(respb4),]

res2<-resmb[complete.cases(resmb),]
res3<-respb[complete.cases(respb),]
#res3<-respb5[complete.cases(respb5),]

common_mp <- data.frame()
common_mp <- as.data.frame(merge(res2, res3, by.x='Human.gene.name', by.y='Human.gene.name'))
pthr=0.1; thr=1
mp<-get_stats(common_mp, "mouse", "pig", pthr, 0.5, thr, thr, 30)
mp<-get_stats(common_mp, "mouse", "pig", pthr, 0.1, thr, thr, 0)
mp[[1]]
cons_mpb<-mp[[2]]


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


