#------------------------------------------------
# Calling lineage-specific methylation
#------------------------------------------------

# The input alignments here are the worker alignments to either the paretnal male or queen 
# i.e. m08 and q08 are the same worker sample but aligned to either the father or mother of those pooled individuals
setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/Bumblebee_files/methylkit_inputs")
library(methylKit)
library(readr)

# --- for each cross

sample.list <- list("w08_male_merged_CpG_evidence.cov", "w08_queen_merged_CpG_evidence.cov",
                    "w19_male_merged_CpG_evidence.cov", "w19_queen_merged_CpG_evidence.cov")

CPGRaw <- methRead(sample.list, 
                   sample.id = list("w08_male","w08_queen","w19_male","w19_queen"),
                   assembly="bter_1.0",
                   treatment=c(0,1,1,0),
                   context="CpG",
                   dbtype = NA,
                   pipeline = "bismarkCoverage",
                   header = T, sep = "\t", mincov=1,
                   dbdir= path)


sample.list <- list("w23_male_merged_CpG_evidence.cov", "w23_queen_merged_CpG_evidence.cov",
                    "w37_male_merged_CpG_evidence.cov", "w37_queen_merged_CpG_evidence.cov")

CPGRaw <- methRead(sample.list, 
                   sample.id = list("w23_male", "w23_queen","w37_male","w37_queen"),
                   assembly="bter_1.0",
                   treatment=c(0,1,1,0),
                   context="CpG",
                   dbtype = NA,
                   pipeline = "bismarkCoverage",
                   header = T, sep = "\t", mincov=1,
                   dbdir= path)

# Filter by coverage NOTE: check how much data we lose here, remember the reads are split between two genomes
filtered_data <- filterByCoverage(CPGRaw,lo.count=8,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)

# Select only CpGs found in all alignments
meth_all_data <- unite(filtered_data, destrand=F) 
nrow(meth_all_data) 
# cross 1 (N genome) 6091
# cross 2 (N genome) 4883

## -------------------------------------------------------------------------

# Filter sites using a binomial test so only keep CpGs which are methylated in at least one sample
df_meth_all <- getData(meth_all_data)

a <- df_meth_all[,5:6]
b <- df_meth_all[,8:9]
c <- df_meth_all[,11:12]
d <- df_meth_all[,14:15]
#e <- df_meth_all[,17:18]
#f <- df_meth_all[,20:21]
#g <- df_meth_all[,23:24]
#h <- df_meth_all[,26:27]

# NOTE: p shouold be the average non-conversion rate (proportion of methylated Cs compared to non-meth Cs)
# So if 1000 methylated Cs compared to 200,000 T's then 1000/200,000 = 0.005
# for a paper: 'the success probability is the non-conversion rate'
bt <- function(a, b, p = 0.005) {binom.test(a, b, 0.005, alternative="greater")$p.value}
allrows <- data.frame(CT=numeric(),Ccount=numeric(),pVal=numeric(),FDR=numeric(),row=numeric())

#for (df in list(a,b,c,d,e,f,g,h)) {
  for (df in list(a,b,c,d)) {
  colnames(df) <- c("CT", "Ccount")
  df$pVal <- mapply(bt, df$Ccount, df$CT)
  df$FDR <- p.adjust(df$pVal, method = "BH", n = length(df$pVal))
  df$row <- as.numeric(rownames(df))
  dfmeth <- subset(df, df$FDR < 0.05)
  allrows <- rbind(allrows, dfmeth)
}
meth_positions <- as.vector(as.numeric(unique(allrows$row))) 
length(meth_positions) 
# cross 2 (N genome): 0 increased to 3
# cross 1 (N genome): 0 increased to 2


subset_methBase <- methylKit::select(meth_all_data, meth_positions)
head(subset_methBase)

# Diff meth between worker alignments to parenta,l genomes
covariates <- data.frame(colony=c("male","female","male","female"))

diff_meth <- calculateDiffMeth(subset_methBase, covariates=covariates, mc.cores = 1)
write.csv(diff_meth, file="all_results_diff_meth.csv")

diff_meth_5 <- getMethylDiff(diff_meth, difference=5, qvalue=0.05)
write.csv(diff_meth_5, file="DMRs_min5percentDiff_qval0.05_MSCfilter.csv")

# cross 1 0/3 significant
# cross 2 0/2 significant