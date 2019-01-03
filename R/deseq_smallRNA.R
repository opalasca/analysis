setwd("~/Desktop/IBD/BGI_analysis")

source("R/deseq_functions.R")
#theme_set(theme_minimal() )  # pre-set the bw theme.

m_new_mirna_file="data/new_miRNA_mus_musculus.csv"
p_new_mirna_file="data/new_miRNA_sus_scrofa.csv"
m_mirna_seq="data/mir_name_sequence_mus_musculus.csv"
p_mirna_seq="data/mir_name_sequence_sus_scrofa.csv"

m_new <- read.csv(m_new_mirna_file, header=TRUE, sep="\t")
p_new <- read.csv(p_new_mirna_file, header=TRUE, sep="\t")
m_seq <- read.csv(m_mirna_seq, header=FALSE, sep="\t")
p_seq <- read.csv(p_mirna_seq, header=FALSE, sep="\t")
colnames(m_seq)=c("id","seq")
colnames(p_seq)=c("id","seq")

m_selected <- filter_overlapping_mirs(m_seq)
p_selected <- filter_overlapping_mirs(p_seq)

miRNAs=c("miRNA")

# Mouse colon #
###############
#data<-get_data_miRNA_v2("data/mir_g1/miRNAs_expressed_all_samples_1545223185.csv", "data/config_mouse.txt", 24, m_selected)
data<-get_data_miRNA_v2("data/mir_g1_s_option/miRNAs_expressed_all_samples_1546186982.csv", "data/config_mouse.txt", 24, m_selected)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
ddsm <- get_deseq(cts2, coldata, colon_samples, 1, 0)
#vsd <- varianceStabilizingTransformation(ddsm, blind=TRUE)
rld <- rlog(ddsm, blind=TRUE)
meanSdPlot(assay(rld))
basic_plots(ddsm, rld, "mouse", "colon","small")
resmcs<-NULL
resmcs <- as.data.frame(get_results(ddsm, "mouse", "colon", "small", 1, mh_orth, pm_orth, m_biotypes, m_seq, "", NULL, NULL))
ddsmcs<-ddsm
hist(resmcs$pvalue)
sig<-return_sig(resmcs, NULL, NULL, 0, 0.05, 1, miRNAs)
clust<-heatmap_DE(sig,rld,"mouse","colon","small","miRNAs","")



# Mouse blood, all samples (timepoint and condition combined and all controls considered as DSS day 0
data<-get_data_miRNA_v2("data/mir_g1_s_option/miRNAs_expressed_all_samples_1546186982.csv", "data/config_mouse.txt", 24, m_selected)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
blood_samples=rownames(coldata[coldata$tissue=="blood",])
ddsm <- get_deseq_time_cond_merged(cts, coldata, blood_samples, 1, 0,"day0_DSS")
rld <- rlog(ddsm, blind=TRUE)
meanSdPlot(assay(rld))
basic_plots(ddsm, rld, "mouse", "blood", "small")
resmbs<-NULL
resmbs <- as.data.frame(get_results(ddsm, "mouse", "blood", "small", 1, mh_orth, pm_orth, 
                                    m_biotypes, m_seq,"","day8_DSS","day0_DSS"))
sig<-return_sig(resmbs, NULL, NULL, 0, 0.05, 1, miRNAs)
heatmap_DE(sig,rld,"mouse","blood","small","miRNAs","t8")
hist(resmbs$pvalue,breaks = 0:20/20)
boxplot(assay(rld))

clust<-heatmap_DE_ts(sig,rld,"mouse","blood","small","miRNAs","t8")
hist(resmbs2$pvalue)
clust<-as.data.frame(clust)
subclust <- clust[clust$cluster=="5",]
for (i in rownames(subclust)){
  print(i)
}
#print(rownames(clust[clust$cluster=="5",]))
kable(rownames(clust[clust$cluster=="5",]))
#plotCounts(ddsmbs, gene="mmu-mir-345-5p", intgroup="condition", returnData=TRUE)



# Pig colon
#data<-get_data_miRNA_v2("data/mir_W_g0/miRNAs_expressed_all_samples_1544544304.csv", "data/config_pig.txt", 30 , p_selected)
data<-get_data_miRNA_v2("data/mir_g1_s_option/miRNAs_expressed_all_samples_1546252867.csv", "data/config_pig.txt", 30 , p_selected)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
#colon_samples=rownames(coldata[coldata$tissue=="colon" & coldata$subject!="P20" & coldata$subject!="P17",]) 
ddsp <- get_deseq_batch(cts, coldata, colon_samples, 1, 0)
#vsd <- varianceStabilizingTransformation(ddsp, blind=TRUE)
rld <- rlog(ddsp, blind=TRUE)
meanSdPlot(assay(rld))
basic_plots(ddsp, rld, "pig", "colon","small")
#resOrdered_uc_co <- get_results(dds, "pig", "colon", "small", 0.1)
respcs<-NULL
respcs <- as.data.frame(get_results(ddsp, "pig", "colon", "small", 1, ph_orth, pm_orth, p_biotypes,p_seq, "",NULL,NULL))
ddspcs<-ddsp
res<-respcs
sig<-return_sig(respcs, NULL, NULL, 0, 0.05, 1, miRNAs)
heatmap_DE(sig,rld,"pig","colon","small","miRNAs","")
hist(respcs$pvalue,breaks = 0:20/20)

# Pig blood all samples (combine timepoint and condition) 
data<-get_data_miRNA_v2("data/mir_g1_s_option/miRNAs_expressed_all_samples_1546252867.csv", "data/config_pig.txt", 30 , p_selected)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$subject!="P17" & coldata$subject!="P15",])
ddsp <- get_deseq_time_cond_merged_plus_batch(cts, coldata, blood_samples, 1, 0, "day0_DSS")
rld <- rlog(ddsp, blind=TRUE)
meanSdPlot(assay(rld))
basic_plots(ddsp, rld, "pig", "blood", "small")
respbs4<-NULL
respbs4 <- as.data.frame(get_results(ddsp, "pig", "blood", "small", 1, ph_orth, pm_orth, 
                                    p_biotypes, p_seq,"","day4_DSS","day0_DSS"))
#respbs2 <- as.data.frame(get_results(ddsp, "pig", "blood", "small", 1, ph_orth, pm_orth, 
 #                                   p_biotypes, p_seq,"","day2_DSS","day0_DSS"))
# P value correction
respbs4_corr <- get_results(ddsp, "pig", "blood", "small", 1, ph_orth, pm_orth, 
                                   p_biotypes, p_seq,"","day4_DSS","day0_DSS")
DESeq2Res<-respbs4_corr
DESeq2Res <- DESeq2Res[, -which(names(DESeq2Res) == "padj")]
FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = F)
FDR.DESeq2Res$param[1, "sd"]
DESeq2Res[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
DESeq2Res[,"pvalue"]  <- FDR.DESeq2Res$pval
hist(FDR.DESeq2Res$pval,breaks = 0:20/20)
hist(DESeq2Res$pval,breaks = 0:20/20)
respbs4_corr<- as.data.frame(DESeq2Res)

sig<-return_sig(respbs4_corr,NULL,  NULL, 0, 0.05,1, miRNAs)
heatmap_DE(sig,rld,"pig","blood","small","miRNAs","t4")
hist(respbs$pvalue,breaks = 0:20/20)




# Mouse blood time series DSS treated 
#data<-get_data_miRNA_v2("data/mir_g1/miRNAs_expressed_all_samples_1545223185.csv", "data/config_mouse.txt", 24, m_selected)
data<-get_data_miRNA_v2("data/mir_g1_s_option/miRNAs_expressed_all_samples_1546186982.csv", "data/config_mouse.txt", 24, m_selected)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" ,])
#blood_samples=rownames(coldata[coldata$tissue=="blood",])
#ddsm <- get_deseq_blood_time(cts, coldata, blood_samples, 1)
#ddsm <- get_deseq_blood_time_full(cts, coldata, blood_samples, 1)
ddsm <- get_deseq_blood_time_no_subject(cts, coldata, blood_samples, 1)
vsd <- varianceStabilizingTransformation(ddsm, blind=TRUE)
rld <- rlog(ddsm, blind=TRUE)
meanSdPlot(assay(vsd))
basic_plots(ddsm, vsd, "mouse", "blood", "total",TRUE)
res1 <- results(ddsm, independentFiltering = TRUE)
resultsNames(ddsm)
res1 <- as.data.frame(results(ddsm, name="time_8_vs_0", test="Wald"))
res1 <- as.data.frame(results(ddsm, independentFiltering = TRUE, contrast=c("time","8","0")))
names(res1)[names(res1) == "log2FoldChange"] = "logFC"
names(res1)[names(res1) == "baseMean"] = "avg.expr"
sig<-return_sig(res1,NULL,NULL, 0, 0.05, 1,miRNAs)
hist(res1$pvalue)
#res2 <- as.data.frame(results(ddsm, independentFiltering = TRUE, contrast=c("time","2","0")))
resmbs8 <- get_results(ddsm, "mouse", "blood","small", 1, mh_orth, mp_orth, m_biotypes,m_seq,"8")
resmbs2 <- get_results(ddsm, "mouse", "blood","small", 1, mh_orth, mp_orth, m_biotypes,m_seq,"2")
sig<-return_sig(resmbs8,resmbs2,NULL, 0, 0.05, 1,miRNAs)
heatmap_DE_ts(sig,rld,"mouse","blood","small","miRNAs")
clust<-heatmap_DE_ts(sig,rld,"mouse","blood","small","miRNAs")
hist(resmbs2$pvalue)
clust<-as.data.frame(clust)
subclust <- clust[clust$cluster=="5",]
for (i in rownames(subclust)){
  print(i)
}
#print(rownames(clust[clust$cluster=="5",]))
kable(rownames(clust[clust$cluster=="5",]))

# Mouse blood DSS vs Control day 8
data<-get_data_miRNA_v2("data/mir_g1_s_option/miRNAs_expressed_all_samples_1546186982.csv", "data/config_mouse.txt", 24, m_selected)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
blood_samples=rownames(coldata[(coldata$time=="8") & coldata$tissue=="blood",])
ddsm <- get_deseq(cts, coldata, blood_samples, 1, 0)
#vsd <- varianceStabilizingTransformation(ddsm, blind=FALSE)
rld <- rlog(ddsm, blind=TRUE)
meanSdPlot(assay(rld))
basic_plots(ddsm, rld, "mouse", "blood", "small")
resmbs<-NULL
resmbs <- as.data.frame(get_results(ddsm, "mouse", "blood", "small", 1, mh_orth, pm_orth, m_biotypes, m_seq, ""))
ddsmbs<-ddsm
res<-resmbs
#lfcthr=0; pthr=0.05
#sig <- res[abs(res$logFC)>lfcthr & res$pvalue<pthr, ]
sig<-return_sig(resmbs, NULL, NULL, 0, 0.05, 1,miRNAs)
sig<-return_sig(resmbs, NULL, NULL, 0, 1, 1,miRNAs)
heatmap_DE(sig,rld,"mouse","blood","small","miRNAs","t8")
hist(resmbs$pvalue)

# Pig blood DSS vs Ctrl, day 4
data<-get_data_miRNA_v2("data/mir_g1_s_option/miRNAs_expressed_all_samples_1546252867.csv", "data/config_pig.txt", 30 , p_selected)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
blood_samples=rownames(coldata[(coldata$time=="4") & coldata$tissue=="blood",])
ddsp <- get_deseq(cts, coldata, blood_samples, 1, 0)
vsd <- varianceStabilizingTransformation(ddsp, blind=FALSE)
rld <- rlog(ddsp, blind=TRUE)
meanSdPlot(assay(vsd))
basic_plots(ddsp, vsd, "mouse", "blood", "small")
respbs<-NULL
respbs <- as.data.frame(get_results(ddsp, "pig", "blood", "small", 1, mh_orth, pm_orth, m_biotypes, m_seq, "",NULL,NULL))
ddspbs<-ddsp
res<-respbs
#lfcthr=0; pthr=0.05
#sig <- res[abs(res$logFC)>lfcthr & res$pvalue<pthr, ]
sig<-return_sig(respbs, NULL, NULL, 0, 0.05, 1,miRNAs)
heatmap_DE(sig,rld,"pig","blood","small","miRNAs")
hist(respbs$pvalue,breaks = 0:20/20)
hist(respbs$pvalue[respbs$avg.expr > 1], breaks = 0:20/20, col = "grey50", border = "white")
#data<-get_data_miRNA("data/miRNAs_expressed_all_samples_1542279087.csv", "data/config_pig.txt", 30)
#data<-get_data_miRNA("data/miRNAs_expressed_all_samples_05_07_2018_t_13_32_49.csv", "data/config_pig.txt", 30)

# Pig blood DSS 
data<-get_data_miRNA_v2("data/mir_g1_s_option/miRNAs_expressed_all_samples_1546252867.csv", "data/config_pig.txt", 30 , p_selected)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS",])
#blood_samples=rownames(coldata[(coldata$time=="4") & coldata$tissue=="blood",])
#ddsp <- get_deseq(cts, coldata, blood_samples, 3,0)
ddsp <- get_deseq_blood_time(cts, coldata, blood_samples, 1 )
vsd <- varianceStabilizingTransformation(ddsp, blind=TRUE)
rld <- rlog(ddsp, blind=TRUE)
meanSdPlot(assay(vsd))
basic_plots(ddsp, vsd, "pig", "blood","small")
respbs4 <- as.data.frame(get_results(ddsp, "pig", "blood", "small", 1, ph_orth, pm_orth, p_biotypes,p_seq,"4",NULL,NULL))
respbs5 <- as.data.frame(get_results(ddsp, "pig", "blood", "small", 1, ph_orth, pm_orth, p_biotypes,p_seq,"5",NULL,NULL))
respbs2 <- as.data.frame(get_results(ddsp, "pig", "blood", "small", 1, ph_orth, pm_orth, p_biotypes,p_seq,"2",NULL,NULL))
hist(respbs4$pvalue)
sig<-return_sig(respbs4,respbs2,respbs5, 0, 0.05, 1,miRNAs)
heatmap_DE_ts(sig,rld,"pig","blood","small","miRNAs")


