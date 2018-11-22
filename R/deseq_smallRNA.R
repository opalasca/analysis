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

colnames(m_seq)=c("fullid","seq")
m_seq$id <- tolower(gsub("_.*","",m_seq$fullid))

#data<-get_data_miRNA("data/miRNAs_expressed_all_samples_1536160042.csv", "data/config_mouse.txt", 24)
data<-get_data_miRNA("data/miRNAs_expressed_all_samples_02_07_2018_t_15_45_27.csv", "data/config_mouse.txt", 24)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
ddsm <- get_deseq(cts, coldata, colon_samples, 3, 0)
vsd <- varianceStabilizingTransformation(ddsm, blind=FALSE)
rld <- rlog(ddsm, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(ddsm, vsd, "mouse", "colon","small")
resmcs<-NULL
resmcs <- as.data.frame(get_results(ddsm, "mouse", "colon", "small", 1, mh_orth, pm_orth, m_biotypes, m_seq, ""))
ddsmcs<-ddsm
res<-resmcs
lfcthr=1; pthr=0.05
sig <- res[abs(res$logFC)>lfcthr & res$pvalue<pthr & res$avg.expr > 10,]
heatmap_DE(sig,rld,"mouse","colon","small")
#heatmap_DE(siglinc,rld,"mouse","colon","total_lincRNA")


#data<-get_data_miRNA("data/miRNAs_expressed_all_samples_1536160042.csv", "data/config_mouse.txt", 24)
data<-get_data_miRNA("data/miRNAs_expressed_all_samples_02_07_2018_t_15_45_27.csv", "data/config_mouse.txt", 24)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
blood_samples=rownames(coldata[(coldata$time=="8") & coldata$tissue=="blood",])
ddsm <- get_deseq(cts, coldata, blood_samples, 3, 0)
vsd <- varianceStabilizingTransformation(ddsm, blind=FALSE)
rld <- rlog(ddsm, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(ddsm, vsd, "mouse", "blood", "small")
resmbs<-NULL
resmbs <- as.data.frame(get_results(ddsm, "mouse", "blood", "small", 1, mh_orth, pm_orth, m_biotypes, m_seq, ""))
ddsmbs<-ddsm
res<-resmbs
lfcthr=1; pthr=0.05
sig <- res[abs(res$logFC)>lfcthr & res$pvalue<pthr & res$avg.expr>10,]
heatmap_DE(sig,rld,"mouse","blood","small")

plotCounts(ddsmbs, gene="mmu-mir-345-5p", intgroup="condition", returnData=TRUE)


# Mouse blood time series
data<-get_data_miRNA("data/miRNAs_expressed_all_samples_02_07_2018_t_15_45_27.csv", "data/config_mouse.txt", 24)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="D" ,])
ddsm <- get_deseq_blood_time(cts, coldata, blood_samples, 3)
vsd <- varianceStabilizingTransformation(ddsm, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(ddsm, vsd, "mouse", "blood", "total",TRUE)
resmbs8 <- get_results(ddsm, "mouse", "blood","total", 1, mh_orth, mp_orth, m_biotypes,m_seq,"8")
resmbs2 <- get_results(ddsm, "mouse", "blood","total", 1, mh_orth, mp_orth, m_biotypes,m_seq,"2")


data<-get_data_miRNA("data/miRNAs_expressed_all_samples_1542279087.csv", "data/config_pig.txt", 30)
data<-get_data_miRNA("data/miRNAs_expressed_all_samples_05_07_2018_t_13_32_49.csv", "data/config_pig.txt", 30)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
#colon_samples=rownames(coldata[coldata$tissue=="colon" & coldata$subject!="P20" & coldata$subject!="P17",]) 
ddsp <- get_deseq(cts, coldata, colon_samples, 3, 0)
vsd <- varianceStabilizingTransformation(ddsp, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(ddsp, vsd, "pig", "colon","small")
#resOrdered_uc_co <- get_results(dds, "pig", "colon", "small", 0.1)
respcs<-NULL
respcs <- as.data.frame(get_results(ddsp, "pig", "colon", "small", 1, ph_orth, pm_orth, p_biotypes,p_seq, ""))
ddspcs<-ddsp

data<-get_data_miRNA("data/miRNAs_expressed_all_samples_1542279087.csv", "data/config_pig.txt", 30)
data<-get_data_miRNA("data/miRNAs_expressed_all_samples_05_07_2018_t_13_32_49.csv", "data/config_pig.txt", 30)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
#blood_samples=rownames(coldata[coldata$time=="4" & coldata$tissue=="blood",])
blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="D",])
blood_samples=rownames(coldata[(coldata$time=="4") & coldata$tissue=="blood",])
#ddsp <- get_deseq(cts, coldata, blood_samples, 3,0)
ddsp <- get_deseq_blood_time(cts, coldata, blood_samples, 3 )
vsd <- varianceStabilizingTransformation(ddsp, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(ddsp, vsd, "pig", "blood","small")
respbs4 <- as.data.frame(get_results(ddsp, "pig", "blood", "small", 1, ph_orth, pm_orth, p_biotypes,"4"))
respbs5 <- as.data.frame(get_results(ddsp, "pig", "blood", "small", 1, ph_orth, pm_orth, p_biotypes,"5"))





