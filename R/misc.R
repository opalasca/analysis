
data<-get_data("data/pig_strtie_quant_tpm.csv", "data/config_pig.txt")
tpmp<-data[[1]]; coldatap<-data[[2]]
data<-get_data("data/mouse_strtie_quant_tpm.csv", "data/config_mouse.txt")
tpmm<-data[[1]]; coldatam<-data[[2]]

tpm<-tpmm
org="mouse"
tpm<-tpmp
org="pig"
tpm<-tpm[rowSums( tpm != 0 ) >= 3,] 
tpm <- log(tpm+1)
write.table(tpm, file=paste("results/", org, "_log_TPM.tsv", sep=""), sep=",",row.names=T, col.names=NA,quote = FALSE)
