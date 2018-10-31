
setwd("~/Desktop/IBD/BGI_analysis")

source("R/deseq_functions.R")
#theme_set(theme_minimal() )  # pre-set the bw theme.

m_orth_file="data/mouse_pig_human_orth.tsv"
p_orth_file="data/pig_human_mouse_orth.tsv"
m_orth <- read.csv(m_orth_file, header=TRUE, sep="\t")
p_orth <- read.csv(p_orth_file, header=TRUE, sep="\t")
mh_orth <- get_orth_table(m_orth, 5, c(1,3,4)) 
mp_orth <- get_orth_table(m_orth, 8, c(1,6,7))
ph_orth <- get_orth_table(p_orth, 5, c(1,3,4)) 
pm_orth <- get_orth_table(p_orth, 9, c(1,7,8))

p_biotypes <- read.csv("data/biotypes/pig_biotypes.txt", header=FALSE, sep="\t")
colnames(p_biotypes) <- c("id","biotype")
m_biotypes <- read.csv("data/biotypes/mouse_biotypes.txt", header=FALSE, sep="\t")
colnames(m_biotypes) <- c("id","biotype")

data<-get_data("data/gene_count_matrix_mouse_annot.csv", "data/config_mouse.txt")
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
colon_samples=rownames(coldata[coldata$tissue=="colon",])
#ddsm <- get_deseq_ref_D(cts, coldata, colon_samples, 3)
ddsm<-NULL
ddsm <- get_deseq(cts, coldata, colon_samples, 1)
#vsd <- vst(ddsm, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
vsd <- rlog(ddsm, blind=FALSE);
basic_plots(ddsm, vsd, "mouse", "colon","total")
resmc<-NULL
resmc <- as.data.frame(get_results(ddsm, "mouse", "colon", "total", 1, mh_orth, mp_orth, m_biotypes))

data<-get_data("data/gene_count_matrix_pig_annot.csv", "data/config_pig.txt")
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
#colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
colon_samples=rownames(coldata[coldata$tissue=="colon" & coldata$subject!="P10" & coldata$subject!="P20",]) 
colon_samples=rownames(coldata[coldata$tissue=="colon" & coldata$subject!="P20",]) 
#ddsp <- get_deseq_ref_D(cts, coldata, colon_samples, 3)
ddsp<-NULL
ddsp <- get_deseq_batch(cts, coldata, colon_samples, 1)
vsd <- vst(ddsp, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
vsd <- rlog(ddsp, blind=FALSE);
meanSdPlot(assay(vsd))
basic_plots(ddsp, vsd, "pig", "colon","total")
respc<-NULL
respc <- as.data.frame(get_results(ddsp, "pig", "colon", "total", 1, ph_orth, pm_orth, p_biotypes))


data<-get_data("data/gene_count_matrix_mouse_annot.csv", "data/config_mouse.txt")
cts<-data[[1]]; coldata<-data[[2]]
#blood_samples=rownames(coldata[(coldata$condition=="C" | coldata$time=="3") & coldata$tissue=="blood",])
blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="D" & (coldata$time=="0" | coldata$time=="8"),])
#blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="D" ,])
ddsm <- get_deseq_blood_paired(cts, coldata, blood_samples, 1)
vsd <- vst(ddsm, blind=FALSE) ##rld <- rlog(dds, blind=FALSE); ntd <- normTransform(dds)
meanSdPlot(assay(vsd))
basic_plots(ddsm, vsd, "mouse", "blood", "total")
resmb <- get_results(ddsm, "mouse", "blood","total", 1, mh_orth, mp_orth, m_biotypes)
ddsmb <- ddsm

data<-get_data("data/gene_count_matrix_pig_annot.csv", "data/config_pig.txt")
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
#blood_samples=rownames(coldata[coldata$time=="5" & coldata$tissue=="blood",]) #use this if comparing with colon
blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="D" & 
                                 (coldata$time=="0" | coldata$time=="4"),])
#blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="D",])
ddsp <- NULL
ddsp <- get_deseq_blood_paired(cts, coldata, blood_samples, 1)
vsd <- vst(ddsp, blind=FALSE) 
meanSdPlot(assay(vsd))
basic_plots(ddsp, vsd, "pig", "blood", "total")
respb<-NULL
respb <- as.data.frame(get_results(ddsp, "pig", "blood", "total", 1, ph_orth, pm_orth, p_biotypes))
ddspb <- ddsp







