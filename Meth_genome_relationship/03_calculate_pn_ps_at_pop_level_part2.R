## -------------------------------------------------------------------------
# Calculate the number of syn and non-syn SNPs that actually occur in each sample
## ------------------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/labelled_snps")

library(sqldf)
library(readr)
library(foreach)
library(doParallel)
library(doBy)
library(dplyr)
library(data.table)
library(reshape2)

# Read in sample snp annotated files
snps <- read_delim("bumblebee_pop_snps_with_info.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)
# Gene annotation file
genes_with_start_and_end <- read_delim("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/genes_with_start_and_end.txt", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)


# Calculate weighted meth for each gene for each sample
colnames(snps)[5]<-"info"
snps <- snps[!is.na(snps$info),]

output <- sqldf("SELECT sample.CHROM,
                    sample.POS,
                    sample.REF,
                    sample.ALT,
                    sample.info,
                    annot.chr,
                    annot.start,
                    annot.end,
                    annot.gene_id
                    FROM snps AS sample
                    LEFT JOIN genes_with_start_and_end AS annot
                    ON sample.CHROM = annot.chr
                    AND (sample.POS >= annot.start AND sample.POS<= annot.end)")
  
output <- output[!is.na(output$gene_id),]
output <- output[,-c(1,2,3,4,7,8)]
output$info[output$info %like% "missense"] <- "non_syn"
output$info[output$info %like% "synonymous"] <- "syn"
output <- output[output$info=="non_syn" | output$info =="syn",]
output <- output %>% group_by_all() %>% count
output <- dcast(output, chr + gene_id ~ info)
output[is.na(output)] <- 0
  
write.table(output, file="pop_level_number_syn_nonsyn_SNPs.txt", quote=F, sep="\t", row.names=F)
