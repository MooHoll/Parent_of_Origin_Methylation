#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(echo=TRUE)

library(sqldf)
library(readr)
library(doBy)
library(dplyr)
library(tidyr)

# To run script: Rscript <name_of_script.R> <input.tsv>
#for file in $(ls *ASM_regions.txt) do;Rscript script.R ${file}; done

base <- tools::file_path_sans_ext(basename(args[1]))
base <- gsub("_ASM_regions","",base)
print(base)

file1 <- args[1]
sample <- read.delim(file1, sep ="\t", header=T)

# Fix file
sample <- separate(data = sample, col = region, into = c("chr", "extra","start", "end"), sep = "\\.")
sample <- sample[,-2]
sample$chr <- paste0(sample$chr, ".1")

annot <- read.delim("HB_genes_with_start_end.txt", sep ="\t", header=F)
colnames(annot)<-c("chr","start","end","gene_id")

output <- sqldf("SELECT sample.chr,
                    sample.start,
                    sample.end,
                    sample.asm_score,
                    sample.coverage,
                    annot.chr,
                    annot.start,
                    annot.end,
                    annot.gene_id
                    FROM sample AS sample
                    LEFT JOIN annot AS annot
                    ON sample.chr = annot.chr
                    AND (sample.start >= annot.start AND sample.end <= annot.end)")
output <- output[!is.na(output$gene_id),]
output <- output[,-6]
colnames(output)[6]<-"gene_start"
colnames(output)[7]<-"gene_end"

myfile <- file.path("./", paste0(base,"_","asm_in_genes.txt"))
write.table(output, file=myfile, quote=F, sep="\t", row.names=F, col.names = T)
