# Making weighted methylation level for each gene for queens/males/worker alignments

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Methylation/PoO_meth/scatter_graphs/meth_per_base")

library(readr)
library(sqldf)
library(doBy)

# Read in sample methylation count files
file.list = list.files(("./"),pattern="*final.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)

# Read in gene start and end file (which includes total CpGs per gene)
genes_with_total_cpgs <- read_table2("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Methylation/PoO_meth/scatter_graphs/genes_with_total_cpgs.txt")
colnames(genes_with_total_cpgs)[5] <- "cpg_count"

# Change chromosomes numbers for the alignments to alternate reference genomes
#chrom_numbersnew <- read_table2("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Methylation/PoO_meth/chrom_numbersnew.txt", 
#                                col_names = FALSE)
#colnames(chrom_numbersnew)<- c("chr","chr1")
#workers <- samples[9:16]
#for(i in seq_along(workers)){
#  output <- merge(workers[[i]], chrom_numbersnew, by="chr")
#  output <- output[,-c(1,2)]
#  colnames(output)[6]<-"chr"
#  myfile <- file.path("./", paste0(i,"_","updated_final.txt"))
#  write.table(output, file=myfile, quote=F, sep="\t", row.names=F)
#}


# Weighted methylation of a region/gene = 
# sum of CpGs in region x all C reads / sum of CpGs in region x (all C + all T reads)


for(i in seq_along(samples)){
  samples[[i]]$count_C <- round((samples[[i]]$coverage/100)*samples[[i]]$freqC, 0)
  df <- samples[[i]]
  output <- sqldf("SELECT sg.chr,
                  sg.base,
                  sg.coverage,
                  sg.count_C,
                  fg.chr,
                  fg.start,
                  fg.end,
                  fg.geneID,
                  fg.cpg_count
                  FROM df AS sg
                  LEFT JOIN genes_with_total_cpgs AS fg 
                  ON sg.chr = fg.chr
                  AND (sg.base >= fg.start AND sg.base <= fg.end)")
  output <- output[!is.na(output$geneID),]
  output <- output[,-c(1,2,5,6,7)]
  check <- summaryBy(coverage + count_C ~ geneID + cpg_count, data=output, FUN=sum) 
  check$weightedMeth <- (check$cpg_count*check$count_C.sum)/(check$cpg_count*check$coverage.sum)
  myfile <- file.path("./", paste0(i,"_","with_weighted_meth.txt"))
  write.table(check, file=myfile, quote=F, sep="\t", row.names=F)
}




## Getting list of genes methylated in at least one sample 

library(readr)
melted_weighted_meth_by_sample_type <- read_delim("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Methylation/PoO_meth/scatter_graphs/weighted_meth/melted_weighted_meth_by_sample_type.txt", 
                                                  "\t", escape_double = FALSE, trim_ws = TRUE)

melted_weighted_meth_by_sample_type <- melted_weighted_meth_by_sample_type[,1:4]
meth_genes <- melted_weighted_meth_by_sample_type[(melted_weighted_meth_by_sample_type$male > 0 |
                                                     melted_weighted_meth_by_sample_type$queens > 0 |
                                                     melted_weighted_meth_by_sample_type$workers > 0),]

meth_genes <- meth_genes[!is.na(meth_genes$geneID),]
# 3550 genes

genes_only <- as.data.frame(meth_genes$geneID)
write.table(genes_only, file = "methylated_genes_leuven_meth.txt", sep="\t",
            col.names = T, row.names = F, quote = F)

  

