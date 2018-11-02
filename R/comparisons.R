
# LogFC comparisons between mouse and pig 

#Mouse vs pig colon
res2<-resmc[complete.cases(resmc),]
res3<-respc[complete.cases(respc),]
common_mp <- data.frame()
common_mp <- as.data.frame(merge(res2, res3, by.x='Human.gene.name', by.y='Human.gene.name'))
pthr=0.1; thr=1
mp<-get_stats(common_mp, "mouse", "pig", pthr, 0.1, thr, thr, 0)
mp<-get_stats(common_mp, "mouse", "pig", 0.1, 0.1, thr, thr, 0)
mp[[1]]
cons_mpc<-mp[[2]]

#Mouse vs pig blood

res2<-resmb2[complete.cases(resmb2),]
res3<-respb2[complete.cases(respb2),]

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


