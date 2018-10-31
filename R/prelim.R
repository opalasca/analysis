setwd("~/Desktop/IBD/BGI_analysis")
source("R/deseq_functions.R")

library(ggplot2)
library(tidyr)
library(purrr)
library(lattice)
library(wesanderson)

draw_density <- function(counts, from, to, genes, name){
  #genes<-pc_mouse
  #from=1
  #to=10
  counts <- log(counts[,from:to]+1)
  counts$id <- rownames(counts)
  d<-as.data.frame(merge(counts, genes, by.x='id', by.y='V1'))
  d %>%
    keep(is.numeric) %>%                     # Keep only numeric columns
    gather() %>%                             # Convert to key-value pairs
    ggplot(aes(value)) +                     # Plot the values
    facet_wrap(~ key)+ #, scales = "free") +   # In separate panels
    geom_density()   
  ggsave(paste("figures/", "density_plots", name, ".pdf", sep=""))
}

count_reads <- function(counts, from, to, genes){
  #counts <- pcts;genes<-pc_pig; from=1; to=30;
  counts$id <- rownames(counts)
  d<-as.data.frame(merge(counts, genes, by.x='id', by.y='V1'))
  sums <-  colSums(d[,(from+1):(to+1)])
  return(sums)
}

library(RColorBrewer)
get_colors<-function(n){
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  pie(rep(1,n), col=sample(col_vector, n))
  return(col_vector)
}

n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

data<-get_data("data/gene_count_matrix_mouse_annot.csv", "data/config_mouse.txt")
mcts<-as.data.frame(data[[1]]); mcoldata<-data[[2]]
data<-get_data("data/gene_count_matrix_pig_annot.csv", "data/config_pig.txt")
pcts<-as.data.frame(data[[1]]); pcoldata<-data[[2]]
p_biotypes <- read.csv("data/biotypes/pig_biotypes.txt", header=FALSE, sep="\t")
colnames(p_biotypes) <- c("id","biotype")
m_biotypes <- read.csv("data/biotypes/mouse_biotypes.txt", header=FALSE, sep="\t")
colnames(m_biotypes) <- c("id","biotype")

counts <- pcts;genes<-p_biotypes; from=1; to=30;
counts$id <- rownames(counts)
d<-as.data.frame(merge(counts, genes, by.x='id', by.y='id'))
psums <- aggregate(. ~ biotype, data=d[,(from+1):(to+2)], FUN=sum)
#barplot(psums,by=biotype)
rownames(psums)<-psums$biotype
psums <- as.matrix(psums[,-1])
psums = psums[ rowSums(psums)>0, ] 
cols<-get_colors(nrow(psums))

pdf(paste("figures/", "biotypes_pig", ".pdf", sep=""))
barchart(t(psums), xlabel="read counts", main="Reads assigned to different gene biotypes", 
         col=cols,
         auto.key=list(space="top", columns=2, title="biotypes", cex.title=1,cex=0.6),
         par.settings=list(superpose.polygon=list(col=cols))
)
#barchart(t(psums), xlabel="read counts", main="Reads assigned to different gene biotypes", 
#         auto.key=list(space="top", columns=5, title="biotypes", cex.title=1))
dev.off()


counts <- mcts;genes<-m_biotypes; from=1; to=24;
counts$id <- rownames(counts)
d<-as.data.frame(merge(counts, genes, by.x='id', by.y='id'))
msums <- aggregate(. ~ biotype, data=d[,(from+1):(to+2)], FUN=sum)
rownames(msums)<-msums$biotype
msums <- as.matrix(msums[,-1])
msums = msums[ rowSums(msums)>10000, ] 
cols<-get_colors(nrow(msums))

pdf(paste("figures/", "biotypes_mouse", ".pdf", sep=""))
barchart(t(msums), xlabel="read counts", main="Reads assigned to different gene biotypes", 
         col=cols,
         auto.key=list(space="top", columns=2, title="biotypes", cex.title=1,cex=0.6),
         par.settings=list(superpose.polygon=list(col=cols))
         )
dev.off()





pc_mouse <- read.csv('data/biotypes/protein_coding_mouse.txt', sep="\t", header=FALSE)
pc_pig <- read.csv('data/biotypes/protein_coding_pig.txt', sep="\t", header=FALSE)
draw_density(mcts,1,18,pc_mouse,"pc_mouse_blood")
draw_density(mcts,19,24,pc_mouse,"pc_mouse_colon")
draw_density(pcts,1,18,pc_pig,"pc_pig_blood")
draw_density(pcts,25,30,pc_pig,"pc_pig_colon")
mb <- count_reads(mcts,1,18,pc_mouse)
mc <- count_reads(mcts,19,24,pc_mouse)
pb <- count_reads(pcts,1,24,pc_pig)
pc <- count_reads(pcts,25,30,pc_pig)
pdf(paste("figures/", "barplot_reads_protein_coding", ".pdf", sep=""))
reads<-c(mc,mb,pc,pb)
barchart(reads,xlim=c(0,90000000), xlabel="read counts", main="Reads assigned to protein coding genes")
dev.off()

