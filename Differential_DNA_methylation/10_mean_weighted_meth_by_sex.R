## -------------------------------------------------------------------------
# Take average weighted methylation level of feature across bio replicates
## -------------------------------------------------------------------------
setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/Differential_methylation/weighted_meth")
library(sqldf)
library(readr)
library(doBy)
library(dplyr)
library(reshape2)

# Make one file covering all samples
file.list = list.files(("./"),pattern="*features.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)


# Make one dataframe for each population
males <- samples[1:4]
queens <- samples[5:8]
workers <- samples[9:12]

for(i in seq_along(males)){
  males[[i]]$caste <- "male"
}
males_all <- as.data.frame(bind_rows(males))
males_merged <- summaryBy(weightedMeth ~ feature + gene_id + start + end +
                        number + caste, data = males_all, FUN=mean)

for(i in seq_along(queens)){
  queens[[i]]$caste <- "queen"
}
queens_all <- as.data.frame(bind_rows(queens))
queens_merged <- summaryBy(weightedMeth ~ feature + gene_id + start + end +
                            number + caste, data = queens_all, FUN=mean)

for(i in seq_along(workers)){
  workers[[i]]$caste <- "worker"
}
workers_all <- as.data.frame(bind_rows(workers))
workers_merged <- summaryBy(weightedMeth ~ feature + gene_id + start + end +
                            number + caste, data = workers_all, FUN=mean)

all_data <- rbind( males_merged, queens_merged, workers_merged)
all_data2 <- dcast(all_data, feature + gene_id + start + end +number ~ caste, value.var = "weightedMeth.mean")


write.table(all_data2, file="weighted_meth_annotation_by_caste.txt", col.names = T,
            row.names = F, quote = F, sep = "\t")



