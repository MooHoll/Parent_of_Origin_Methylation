# Merging pvalues and frequency information 

setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022")
library(readr)

male_high_p <- read_delim("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/GO_analysis/outputs/M_hyper_vsWQ.txt", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE)

male_high_freq <- read_delim("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/GO_analysis/outputs/REVIGO/M_hyper_vsWQ.tsv", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE)
male_high_freq <- male_high_freq[,c(1,2,5)]
colnames(male_high_freq)[1] <- "GOBPID"

both <- merge(male_high_p, male_high_freq, by = "GOBPID")

write.table(both, file="info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)
