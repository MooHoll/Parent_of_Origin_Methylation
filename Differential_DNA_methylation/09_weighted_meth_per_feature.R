## -------------------------------------------------------------------------
# Weighted methylation per annotation for each sex
## -------------------------------------------------------------------------

library(sqldf)
library(readr)
library(doBy)
library(dplyr)
library(foreach)
library(doParallel)

# Read in sample methylation count files
file.list = list.files(("./"),pattern="*final_coverage.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)
sample_names <- list("m08","m19","m23","m37","q08","q19","q23","q37","w08","w19","w23","w37")
names(samples) <- sample_names

# Read in gene with start/end and total CpGs per gene
annotation_with_total_cpgs <- read_table2("ref_Bter_1.0_top_level_numbered_exons_with_total_cpgs.txt")
colnames(annotation_with_total_cpgs)[7] <- "cpg_count"
colnames(annotation_with_total_cpgs)[2] <- "chr"

registerDoParallel(cores = 4)

# Calculate weighted meth for each gene for each sample
foreach(i = seq_along(samples)) %dopar% {
  df <- samples[[i]]
  df <- subset(df, total_coverage > 5)
  output <- sqldf("SELECT sample.chr,
                    sample.cpg,
                    sample.count_c,
                    sample.total_coverage,
                    annot.chr,
                    annot.start,
                    annot.end,
                    annot.gene_id,
                    annot.number,
                    annot.feature,
                    annot.cpg_count
                    FROM df AS sample
                    LEFT JOIN annotation_with_total_cpgs AS annot
                    ON sample.chr = annot.chr
                    AND (sample.cpg >= annot.start AND sample.cpg <= annot.end)")
  output <- output[!is.na(output$gene_id),]
  output <- output[,-c(1,2)]
  check <- summaryBy(total_coverage + count_c ~ chr + feature + gene_id + start + end + number + cpg_count, data=output, FUN=sum) 
  check$weightedMeth <- (check$cpg_count*check$count_c.sum)/(check$cpg_count*check$total_coverage.sum)
  myfile <- file.path("./", paste0(names(samples[i]),"_","weighted_meth_all_features.txt"))
  write.table(check, file=myfile, quote=F, sep="\t", row.names=F)
}


## -------------------------------------------------------------------------
