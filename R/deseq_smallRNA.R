setwd("~/Desktop/IBD/BGI_analysis")

source("R/deseq_functions.R")
#theme_set(theme_minimal() )  # pre-set the bw theme.

m_new_mirna_file="data/new_miRNA_mus_musculus.csv"
p_new_mirna_file="data/new_miRNA_sus_scrofa.csv"

m_new <- read.csv(m_new_mirna_file, header=TRUE, sep="\t")
p_new <- read.csv(p_new_mirna_file, header=TRUE, sep="\t")
m_seq <- read.csv("data/mir_name_sequence_mus_musculus.csv", header=FALSE, sep="\t")
p_seq <- read.csv("data/mir_name_sequence_sus_scrofa.csv", header=FALSE, sep="\t")
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
coldata<-coldata[order(coldata$condition,coldata$time),]
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
ddsm <- get_deseq(cts, coldata, colon_samples, 1, 0)
#vsd <- varianceStabilizingTransformation(ddsm, blind=TRUE)
rldc <- rlog(ddsm, blind=TRUE)
meanSdPlot(assay(rld))
#rld <- rld[,colnames]
basic_plots(ddsm, rldc, "mouse", "colon","small")
resmcs<-NULL
resmcs <- as.data.frame(get_results(ddsm, "mouse", "colon", "small", 1, mh_orth, pm_orth, m_biotypes,"", m_seq, "", NULL, NULL))
ddsmcs<-ddsm
hist(resmcs$pvalue)
rlog_mouse_colon_small <- assay(rldc)
meanrldcs<-get_mean(rlog_mouse_colon_small, coldata, "mouse", "colon", "small")
draw_volcano("mouse","colon", "small",resmcs, 0.1, 1,"")
sig<-return_sig(resmcs, NULL, NULL, 0.5, 0.05, 0.1, miRNAs)
rld<-rldc; res<-resmcs
rows <- match(sig$id, row.names(rld))
mat <- assay(rld)[rows,]
mat<-t(scale(t(mat)))
#heatmap_DE_kmeans(sig,rld, mat,"mouse","colon","small","miRNAs",3, FALSE)
set.seed(8)
km<- kmeans(mat,3,iter.max=1000) # determine how many clusters you want, I specify 2 here
m.kmeans<- cbind(mat, km$cluster) # combine the cluster with the matrix
o <- order(m.kmeans[,dim(m.kmeans)[2]]) # order the last column
m.kmeans <- m.kmeans[o,] # order the matrix according to the order of the last column
colnames(m.kmeans)[dim(m.kmeans)[2]] <- c("cluster")
heatmap_DE_kmeans(sig,rld,m.kmeans[,1:6],"mouse","colon","small","miRNAs",3, FALSE)
m.kmeans<-as.data.frame(m.kmeans)
kable(rownames(m.kmeans[m.kmeans$cluster=="1",]))
#plotCounts(ddsmbs, gene="mmu-mir-345-5p", intgroup="condition", returnData=TRUE)
#heatmap_DE_rownames(sig,rld,"mouse","colon","small","miRNAs","day8", 2)



# Mouse blood, all samples (timepoint and condition combined and all controls considered as DSS day 0
data<-get_data_miRNA_v2("data/mir_g1_s_option/miRNAs_expressed_all_samples_1546186982.csv", "data/config_mouse_rename.txt", 24, m_selected)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
coldata<-coldata[order(coldata$condition,coldata$time),]
blood_samples=rownames(coldata[coldata$tissue=="blood",])
ddsm <- get_deseq_time_cond_merged(cts, coldata, blood_samples, 1, 0,"day0_DSS")
rldb <- rlog(ddsm, blind=TRUE)
#filtered_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" & coldata$time != "2",])
#filtered_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" ,])
#rld <- rld[,colnames(rld) %in% filtered_samples]
meanSdPlot(assay(rldb))
basic_plots(ddsm, rldb, "mouse", "blood", "small",TRUE)
resmbs<-NULL
resmbs <- as.data.frame(get_results(ddsm, "mouse", "blood", "small", 1, mh_orth, pm_orth, 
                                    m_biotypes,"", m_seq,"","day8_DSS","Control"))
resmbs2 <- as.data.frame(get_results(ddsm, "mouse", "blood", "small", 1, mh_orth, pm_orth, 
                                    m_biotypes,"", m_seq,"","day2_DSS","Control"))
#sig<-return_sig(resmbs, resmbs2, NULL, 0.5, 0.05, 0.1, miRNAs)
sig<-return_sig(resmbs, NULL, NULL, 0, 0.05, 0.3, miRNAs)
#set.seed(9)
draw_volcano("mouse","blood", "small",resmbs, 0.1, 1,"Day8_vs_control")
draw_volcano("mouse","blood", "small",resmbs2, 0.1, 1,"Day2_vs_control")
hist(resmbs$pvalue,breaks = 0:20/20)
boxplot(assay(rldb))
rlog_mouse_blood_small <- assay(rldb)
meanrldbs<-get_mean(rlog_mouse_blood_small, coldata, "mouse", "blood", "small")
#Compute clustering on a subset of samples only (e.g. DSS samples)
filtered_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" ,])
rldbf <- rldb[,colnames(rldb) %in% filtered_samples]
rld<-rldbf; res<-resmbs
rows <- match(sig$id, row.names(rld))
mat <- assay(rld)[rows,]
mat<-t(scale(t(mat)))
heatmap_DE_kmeans(sig,rld, mat,"mouse","blood","small","miRNAs",3, FALSE)
set.seed(8)
km<- kmeans(mat,5,iter.max=1000) # determine how many clusters you want, I specify 2 here
#Show heatmap for all samples, even if clustering was computed on a subset of samples
rld<-rldb
mat <- assay(rld)[rows,]
mat<-t(scale(t(mat)))
m.kmeans<- cbind(mat, km$cluster) # combine the cluster with the matrix
o <- order(m.kmeans[,dim(m.kmeans)[2]]) # order the last column
m.kmeans <- m.kmeans[o,] # order the matrix according to the order of the last column
colnames(m.kmeans)[dim(m.kmeans)[2]] <- c("cluster")
heatmap_DE_kmeans_ts(sig,rld,m.kmeans[,1:18],"mouse","blood","small","miRNAs",3, FALSE)
m.kmeans<-as.data.frame(m.kmeans)
blood_M72_small<- m.kmeans[m.kmeans$cluster=="1",] 
kable(rownames(m.kmeans[m.kmeans$cluster=="1",]))
plotCounts(ddsm, gene="mmu-mir-100-5p", intgroup="condition", returnData=TRUE)

# Pig colon
#data<-get_data_miRNA_v2("data/mir_W_g0/miRNAs_expressed_all_samples_1544544304.csv", "data/config_pig.txt", 30 , p_selected)
data<-get_data_miRNA_v2("data/mir_g1_s_option/miRNAs_expressed_all_samples_1546252867.csv", "data/config_pig.txt", 30 , p_selected)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
coldata<-coldata[order(coldata$condition,coldata$time),]
colon_samples=rownames(coldata[coldata$tissue=="colon",]) 
#colon_samples=rownames(coldata[coldata$tissue=="colon" & coldata$subject!="P20" & coldata$subject!="P17",]) 
ddsp <- get_deseq_batch(cts, coldata, colon_samples, 1, 0)
#vsd <- varianceStabilizingTransformation(ddsp, blind=TRUE)
rld <- rlog(ddsp, blind=TRUE)
meanSdPlot(assay(rld))
#basic_plots(ddsp, rld, "pig", "colon","small")
#resOrdered_uc_co <- get_results(dds, "pig", "colon", "small", 0.1)
respcs<-NULL
respcs <- as.data.frame(get_results(ddsp, "pig", "colon", "small", 1, ph_orth, pm_orth, p_biotypes,"",p_seq, "",NULL,NULL))
ddspcs<-ddsp
respcscorr <- get_results_corrected_pval(ddsp, "pig", "colon", "small", 1, ph_orth, pm_orth, 
                                          p_biotypes, p_seq,"",NULL,NULL)
sig<-return_sig(respcscorr, NULL, NULL, 0, 0.05, 0.2, miRNAs)
heatmap_DE_rownames(sig,rld,"pig","colon","small","miRNAs","")
hist(respcscorr$pvalue,breaks = 0:20/20)
rlog_pig_colon_small <- assay(rld)
meanrldcsp<-get_mean(rlog_pig_colon_small, coldata, "pig", "colon", "small")


# Pig blood all samples (combine timepoint and condition) 
data<-get_data_miRNA_v2("data/mir_g1_s_option/miRNAs_expressed_all_samples_1546252867.csv", "data/config_pig.txt", 30 , p_selected)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
coldata<-coldata[order(coldata$condition,coldata$time),]
#blood_samples=rownames(coldata[coldata$tissue=="blood" & coldata$subject!="P17" & coldata$subject!="P15",])
blood_samples=rownames(coldata[coldata$tissue=="blood",])
ddsp <- get_deseq_time_cond_merged_plus_batch(cts, coldata, blood_samples, 1, 0, "day0_DSS")
rld <- rlog(ddsp, blind=TRUE)
meanSdPlot(assay(rld))
filtered_samples=rownames(coldata[coldata$tissue=="blood" & coldata$condition=="DSS" ,])
rld <- rld[,colnames(rld) %in% filtered_samples]
#basic_plots(ddsp, rld, "pig", "blood", "small")
respbs4<-NULL
respbs4 <- as.data.frame(get_results(ddsp, "pig", "blood", "small", 1, ph_orth, pm_orth, 
                                    p_biotypes,"", p_seq,"","day4_DSS","day0_DSS"))
#respbs2 <- as.data.frame(get_results(ddsp, "pig", "blood", "small", 1, ph_orth, pm_orth, 
 #                                   p_biotypes, p_seq,"","day2_DSS","day0_DSS"))
# P value correction
respbs4corr <- get_results_corrected_pval(ddsp, "pig", "blood", "small", 1, ph_orth, pm_orth, 
                                   p_biotypes, p_seq,"","day4_DSS","day0_DSS")
respbs2corr <- get_results_corrected_pval(ddsp, "pig", "blood", "small", 1, ph_orth, pm_orth, 
                                   p_biotypes, p_seq,"","day2_DSS","day0_DSS")
sig<-return_sig(respbs4corr,NULL,NULL, 0, 0.05, 0.2, miRNAs)
#sig<-return_sig(respbs4,NULL,NULL, 0, 0.05,0.3, miRNAs)
heatmap_DE_ts(sig,rld,"pig","blood","small","miRNAs")
hist(respbs2corr$pvalue,breaks = 0:20/20)
rlog_pig_blood_small <- assay(rld)
meanrldbsp<-get_mean(rlog_pig_blood_small, coldata, "pig", "blood", "small")



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
data<-get_data_miRNA_v2("data/mir_g1_s_option/miRNAs_expressed_all_samples_1546186982.csv", "data/config_mouse_rename.txt", 24, m_selected)
cts<-NULL; coldata<-NULL
cts<-data[[1]]; coldata<-data[[2]]
blood_samples=rownames(coldata[(coldata$time=="8") & coldata$tissue=="blood",])
ddsm8 <- get_deseq(cts, coldata, blood_samples, 2, 0)
#vsd <- varianceStabilizingTransformation(ddsm, blind=FALSE)
rld <- rlog(ddsm8, blind=TRUE)
meanSdPlot(assay(rld))
basic_plots(ddsm8, rld, "mouse", "blood", "small")
resmbs8 <- as.data.frame(get_results(ddsm8, "mouse", "blood", "small", 1, mh_orth, pm_orth,m_biotypes,"",m_seq, "",NULL,NULL))
ddsmbs8<-ddsm8
sig<-return_sig(resmbs8, NULL, NULL, 0, 0.05, 0.3, miRNAs)
res<-resmbs; 
rows <- match(sig$id, row.names(rld))
mat <- assay(rld)[rows,]
mat<-t(scale(t(mat)))
#heatmap_DE_kmeans(sig,rld, mat,"mouse","blood","small","miRNAs",3, FALSE)
set.seed(8)
km<- kmeans(mat,2,iter.max=1000) # determine how many clusters you want, I specify 2 here
m.kmeans<- cbind(mat, km$cluster) # combine the cluster with the matrix
o <- order(m.kmeans[,dim(m.kmeans)[2]]) # order the last column
m.kmeans <- m.kmeans[o,] # order the matrix according to the order of the last column
colnames(m.kmeans)[dim(m.kmeans)[2]] <- c("cluster")
heatmap_DE_kmeans(sig,rld,m.kmeans[,1:6],"mouse","blood","small","miRNAs",3, FALSE)
hist(resmbs8$pvalue)

sig8<-return_sig(resmbs8, NULL, NULL, 1, 0.05, 0.9, miRNAs)
siga<-return_sig(resmbs, NULL, NULL, 1, 0.05, 0.9, miRNAs)
common_only <- merge(siga, sig8, by="id")
common <- merge(siga, sig8, by="id", all.x=TRUE, all.y=TRUE)
m<- common[,c(1,2,7,3,11,16,12)]

plot_gene_time_mir("mmu-mir-20a-5p",resmbs,ddsm)
plot_gene_time_mir("mmu-let-7i-5p",resmbs,ddsm)


plot_gene("mmu-mir-3963","",ddsm8)

a <- dim(siga)[1]
b <- dim(sig8)[1]
ab <- dim(common_only)[1]
grid.newpage()
plot <- draw.pairwise.venn(a, b, ab, category = c("all samples", "day 8"),cex=2,cat.cex=2, cat.just=list(c(0.2,1), c(1,1)))
pdf("figures/venn_diagram_blood_mirna.pdf",width=6,height=5) 
grid.draw(plot);
dev.off()


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



#Hierarchical clustering code
clust<-heatmap_DE_ts_clust(sig,rld,"mouse","blood","small","miRNAs","t8")
#hist(resmbs2$pvalue)
clust<-as.data.frame(clust)
subclust <- clust[clust$cluster=="3",]
for (i in rownames(subclust)){
  print(i)
}
#print(rownames(clust[clust$cluster=="5",]))
