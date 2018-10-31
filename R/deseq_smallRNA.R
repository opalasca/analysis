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

data<-get_data_miRNA("data/miRNAs_expressed_all_samples_1536160042.csv", "data/config_mouse.txt", 24)
#data<-get_data_miRNA("data/miRNAs_expressed_all_samples_02_07_2018_t_15_45_27.csv", "data/config_mouse.txt", 24)
cts<-data[[1]]; coldata<-data[[2]]
dds <- get_deseq(cts, coldata, c(19:24), 3)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(dds, vsd, "mouse", "colon","small")
resOrdered_uc_co <- get_results(dds, "mouse", "colon", "small", 0.1)

data<-get_data_miRNA("data/miRNAs_expressed_all_samples_1536160042.csv", "data/config_mouse.txt", 24)
#data<-get_data_miRNA("data/miRNAs_expressed_all_samples_02_07_2018_t_15_45_27.csv", "data/config_mouse.txt", 24)
cts<-data[[1]]; coldata<-data[[2]]
blood_samples=rownames(coldata[(coldata$condition=="C" | coldata$time=="3") & coldata$tissue=="blood",])
dds <- get_deseq(cts, coldata, blood_samples, 3)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(dds, vsd, "mouse", "blood", "small")
resOrdered_uc_co <- get_results(dds, "mouse", "blood","small", 0.1)

data<-get_data_miRNA("data/miRNAs_expressed_all_samples_1536167449.csv", "data/config_pig.txt", 30)
cts<-data[[1]]; coldata<-data[[2]]
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
dds <- get_deseq(cts, coldata, colon_samples, 3)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(dds, vsd, "pig", "colon","small")
resOrdered_uc_co <- get_results(dds, "pig", "colon", "small", 0.1)

blood_samples=rownames(coldata[(coldata$condition=="C" | coldata$time=="4") & coldata$tissue=="blood",])
dds <- get_deseq(cts, coldata, blood_samples, 3)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
meanSdPlot(assay(vsd))
basic_plots(dds, vsd, "pig", "blood","small")
resOrdered_uc_co <- get_results(dds, "pig", "blood", "small", 0.1)




