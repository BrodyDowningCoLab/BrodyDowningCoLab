#Andrew Phan
#2/15/21

rm(list=ls()) #clears variables
cat("\014") #clears console

library(data.table)

tab = fread(file = "/home/data/Shared/TCGA/TCGA_COAD_Methylation_Chr21.csv", header = T) #read in methylation data
tab = tab[tab$sample_barcode!="",] #remove rows with blank sample IDs
ids = unique(tab$sample_barcode) #get unique sample IDs


#BED file format:
#Chromosome, chromStart, chromEnd, sample ID, and beta value
makeBED <- function(id) {
  index = tab$sample_barcode==id
  BED = cbind(tab$chromosome[index], tab$position[index], tab$position[index], tab$sample_barcode[index], tab$beta_value[index])
  return(BED)
}

res <- lapply(ids, makeBED)

#write a tab delimited file for each sample ID
toTxt <- function(data){
  outdir = "/home/data/Shared/TCGA/BED_newDownload/"
  sampID = data.frame(data)[1,4] #store sample ID
  fwrite(data, file = paste(outdir, sampID, "_BED.bed", sep = ""), sep = "\t", col.names = F) #bedtools does not accept headers (column names)
  print(paste("Made: ", outdir, sampID, "_BED.bed", sep = ""))
}

for (i in 1:length(res)){
  lapply(res[i], toTxt)
}