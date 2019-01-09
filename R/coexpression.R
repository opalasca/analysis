setwd("~/Desktop/IBD/BGI_analysis")
source("R/deseq_functions.R")
require(gdata)


# Validated targets from mirTarBase
mmu_validated = read.xls("data/mmu_MTI.xls", sheet = 1, header = TRUE)
ssc_validated = read.xls("data/ssc_MTI.xls", sheet = 1, header = TRUE)



# Commands to obtain mirna-host pairs from mirbase download files
mirna_host_mouse = read.csv("data/mirna_context_mouse.txt", header=FALSE, sep="\t")
mirna_host_pig = read.csv("data/mirna_context_pig.txt", header=FALSE, sep="\t")
prec_id_name <- read.csv("data/mirna.txt", header=FALSE, sep="\t")
gene_tr_name_mouse <- read.csv("data/mouse_gene_transcript_name.tsv", header=TRUE, sep="\t")
gene_tr_name_pig <- read.csv("data/mouse_gene_transcript_name.tsv", header=TRUE, sep="\t")
mirna_host <- mirna_host_mouse[,c(1:2,4)]
names(mirna_host)<- c("internal_id","host_tr_id","region")
id_name <- prec_id_name[,c(1,3,4,5)]
names(id_name)<- c("internal_id","name","alternative_name","full_name")
m_host <- as.data.frame(merge(mirna_host, id_name, by="internal_id",all.x=TRUE ))
