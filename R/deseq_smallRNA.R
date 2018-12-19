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
#m_seq$id <- tolower(gsub("_.*","",m_seq$fullid))
colnames(p_seq)=c("id","seq")
#p_seq$id <- tolower(gsub("_.*","",p_seq$fullid))

seqf<-m_seq
seqf<-seqf[order(seqf$seq),]
i=2
while (i<length(seqf$seq)){
  #print(i)
  k=1
  a<-grepl(seqf[i-1,2],seqf[i,2])
  if (a==TRUE){
    b=TRUE;
    k=i;
    while (b==TRUE){
      b<-grepl(seqf[k,2],seqf[k+1,2])
      k=k+1
    }  
    low=i-1;high=k-1;
    if((high-low)>0){
      print(paste(low, high, sep=" "))
      i=k
      print(seqf[low:high,])
      #print(seqf[i-1,])
      #print(seqf[i,])
    } 
  }
  else{
      i=i+1
    } 
}  
  


#data<-get_data_miRNA("data/miRNAs_expressed_all_samples_1536160042.csv", "data/config_mouse.txt", 24)
#data<-get_data_miRNA("data/miRNAs_expressed_all_samples_02_07_2018_t_15_45_27.csv", "data/config_mouse.txt", 24)
#data<-get_data_miRNA_v2("data/mir_W_g0/miRNAs_expressed_all_samples_1544542395.csv", "data/config_mouse.txt", 24)
data<-get_data_miRNA_v2("data/mir_g1/miRNAs_expressed_all_samples_1545223185.csv", "data/config_mouse.txt", 24)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
ddsm <- get_deseq(cts, coldata, colon_samples, 1, 0)
vsd <- varianceStabilizingTransformation(ddsm, blind=TRUE)
rld <- rlog(ddsm, blind=TRUE)
meanSdPlot(assay(vsd))
basic_plots(ddsm, vsd, "mouse", "colon","small")
resmcs<-NULL
resmcs <- as.data.frame(get_results(ddsm, "mouse", "colon", "small", 1, mh_orth, pm_orth, m_biotypes, m_seq, ""))
ddsmcs<-ddsm
res<-resmcs
hist(resmcs$pvalue)

lfcthr=0; pthr=0.05
sig <- res[abs(res$logFC)>lfcthr & res$pvalue<pthr ,]
heatmap_DE(sig,rld,"mouse","colon","small","miRNAs")
#heatmap_DE(siglinc,rld,"mouse","colon","total_lincRNA")


#data<-get_data_miRNA("data/miRNAs_expressed_all_samples_1536160042.csv", "data/config_mouse.txt", 24)
#data<-get_data_miRNA("data/miRNAs_expressed_all_samples_02_07_2018_t_15_45_27.csv", "data/config_mouse.txt", 24)
#data<-get_data_miRNA_v2("data/mir_W_g0/miRNAs_expressed_all_samples_1544542395.csv", "data/config_mouse.txt", 24)
data<-get_data_miRNA_v2("data/mir_g1/miRNAs_expressed_all_samples_1545223185.csv", "data/config_mouse.txt", 24)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
blood_samples=rownames(coldata[(coldata$time=="8") & coldata$tissue=="blood",])
ddsm <- get_deseq(cts, coldata, blood_samples, 1, 0)
vsd <- varianceStabilizingTransformation(ddsm, blind=FALSE)
rld <- rlog(ddsm, blind=TRUE)
meanSdPlot(assay(vsd))
basic_plots(ddsm, vsd, "mouse", "blood", "small")
resmbs<-NULL
resmbs <- as.data.frame(get_results(ddsm, "mouse", "blood", "small", 1, mh_orth, pm_orth, m_biotypes, m_seq, ""))
ddsmbs<-ddsm
res<-resmbs
lfcthr=0; pthr=0.05
sig <- res[abs(res$logFC)>lfcthr & res$pvalue<pthr, ]
heatmap_DE(sig,rld,"mouse","blood","small","miRNAs")
hist(resmcs$pvalue)


#plotCounts(ddsmbs, gene="mmu-mir-345-5p", intgroup="condition", returnData=TRUE)


# Mouse blood time series
#data<-get_data_miRNA("data/miRNAs_expressed_all_samples_02_07_2018_t_15_45_27.csv", "data/config_mouse.txt", 24)
#data<-get_data_miRNA_v2("data/mir_W_g0/miRNAs_expressed_all_samples_1544542395.csv", "data/config_mouse.txt", 24)
data<-get_data_miRNA_v2("data/mir_g1/miRNAs_expressed_all_samples_1545223185.csv", "data/config_mouse.txt", 24)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" ,])
ddsm <- get_deseq_blood_time(cts, coldata, blood_samples, 1)
vsd <- varianceStabilizingTransformation(ddsm, blind=TRUE)
rld <- rlog(ddsm, blind=TRUE)
meanSdPlot(assay(vsd))
basic_plots(ddsm, vsd, "mouse", "blood", "total",TRUE)
resmbs8 <- get_results(ddsm, "mouse", "blood","total", 1, mh_orth, mp_orth, m_biotypes,m_seq,"8")
resmbs2 <- get_results(ddsm, "mouse", "blood","total", 1, mh_orth, mp_orth, m_biotypes,m_seq,"2")
res<-resmbs8
lfcthr=0; pthr=0.05
sig <- res[abs(res$logFC)>lfcthr & res$pvalue<pthr & res$avg.expr>0,]
#heatmap_DE_ts(sig,rld,"mouse","blood","small","miRNAs")
sig<-return_sig(resmbs8,resmbs2,NULL, 0, 0.05)
sig<-return_sig(resmbs2,NULL,NULL, 0, 0.05)

heatmap_DE_ts(sig,rld,"mouse","blood","small","miRNAs")
hist(resmbs2$pvalue)

data<-get_data_miRNA_v2("data/mir_W_g0/miRNAs_expressed_all_samples_1544544304.csv", "data/config_pig.txt", 30)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
#colon_samples=rownames(coldata[coldata$tissue=="colon" & coldata$subject!="P20" & coldata$subject!="P17",]) 
ddsp <- get_deseq(cts, coldata, colon_samples, 1, 0)
vsd <- varianceStabilizingTransformation(ddsp, blind=TRUE)
rld <- rlog(ddsp, blind=TRUE)
meanSdPlot(assay(vsd))
basic_plots(ddsp, vsd, "pig", "colon","small")
#resOrdered_uc_co <- get_results(dds, "pig", "colon", "small", 0.1)
respcs<-NULL
respcs <- as.data.frame(get_results(ddsp, "pig", "colon", "small", 1, ph_orth, pm_orth, p_biotypes,p_seq, ""))
ddspcs<-ddsp
res<-respcs
lfcthr=0; pthr=1; pvthr=0.05
sig <- res[abs(res$logFC)>lfcthr & res$pval<pvthr ,]
heatmap_DE(sig,rld,"pig","colon","small","miRNAs")
#hist(respcs$pvalue)

#data<-get_data_miRNA("data/miRNAs_expressed_all_samples_1542279087.csv", "data/config_pig.txt", 30)
#data<-get_data_miRNA("data/miRNAs_expressed_all_samples_05_07_2018_t_13_32_49.csv", "data/config_pig.txt", 30)
data<-get_data_miRNA_v2("data/mir_W_g0/miRNAs_expressed_all_samples_1544544304.csv", "data/config_pig.txt", 30)
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
respbs4 <- as.data.frame(get_results(ddsp, "pig", "blood", "small", 1, ph_orth, pm_orth, p_biotypes,p_seq,"4"))
respbs5 <- as.data.frame(get_results(ddsp, "pig", "blood", "small", 1, ph_orth, pm_orth, p_biotypes,p_seq,"5"))
respbs2 <- as.data.frame(get_results(ddsp, "pig", "blood", "small", 1, ph_orth, pm_orth, p_biotypes,p_seq,"2"))
hist(respbs4$pvalue)
sig<-return_sig(respbs4,respbs2,respbs5, 0, 0.05)
heatmap_DE_ts(sig,rld,"pig","blood","small","miRNAs")


