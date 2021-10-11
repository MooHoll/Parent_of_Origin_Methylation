#------------------------------------------------
# Calling parent-specific methylation
#------------------------------------------------
# On ALICE2 R/3.6.1
library(methylKit)
library(readr)

## -------------------------------------------------------------------------

# Process bam files to text files for easier downstream messing about
file.list <- list("w08_male_reads.bam","w19_male_reads.bam",
                  "w23_male_reads.bam","w37_male_reads.bam",
                  "w08_queen_reads.bam","w19_queen_reads.bam", 
                  "w23_queen_reads.bam","w37_queen_reads.bam")

raw_data <- processBismarkAln(file.list,
                              sample.id = list("w08_male","w19_male","w23_male","w37_male",
                                               "w08_queen","w19_queen", "w23_queen", "w37_queen"),
                              treatment = c(1,1,1,1,0,0,0,0),
                              assembly="bter_1.0", 
                              read.context="CpG",
                              save.context = "CpG",
                              save.folder=getwd(),
                              save.db = TRUE)

## -------------------------------------------------------------------------

# Read in previously created text files
file.list <- list("w08_male_CpG.txt", "w19_male_CpG.txt", 
                  "w08_queen_CpG.txt", "w19_queen_CpG.txt")

raw_data <- methRead(file.list,
                     sample.id = list("w08_male","w19_male",
                                      "w08_queen","w19_queen"),
                     treatment = c(1,1,0,0),
                     assembly="bter_1.0", 
                     context="CpG")


file.list <- list("w23_male_CpG.txt", "w37_male_CpG.txt", 
                  "w23_queen_CpG.txt", "w37_queen_CpG.txt")

raw_data <- methRead(file.list,
                     sample.id = list("w23_male","w37_male",
                                      "w23_queen","w37_queen"),
                     treatment = c(1,1,0,0),
                     assembly="bter_1.0", 
                     context="CpG")

# For all data together (then only get 3 CpGs shared across all data when unite)
file.list <- list("w08_male_CpG.txt", "w19_male_CpG.txt", "w23_male_CpG.txt", "w37_male_CpG.txt", 
                  "w08_queen_CpG.txt", "w19_queen_CpG.txt", "w23_queen_CpG.txt", "w37_queen_CpG.txt")

raw_data <- methRead(file.list,
                     sample.id = list("w08_male","w19_male","w23_male","w37_male",
                                      "w08_queen","w19_queen", "w23_queen", "w37_queen"),
                     treatment = c(1,1,1,1,0,0,0,0),
                     assembly="bter_1.0", 
                     context="CpG")

## -------------------------------------------------------------------------

# Filter data for outliers and coverage
filtered_data <- filterByCoverage(raw_data,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)

## -------------------------------------------------------------------------

# Only text CpGs present in all samples (maybe too stingent, can adjust later)
meth_all_data <- unite(filtered_data, destrand=TRUE) 
nrow(meth_all_data) 
# first cross 1003
# second cross 557
# both 3

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
# first cross: 1
# second cross: 0
# both: 0
subset_methBase <- methylKit::select(meth_all_data, meth_positions)

## -------------------------------------------------------------------------

# Save the dataframe for later use, including making nice plots
methBase_ob <- getData(subset_methBase)
write.table(methBase_ob, file="MaleReads_vs_QueenReads_objectmethbase.txt", quote=F, row.names = F, sep = '\t')

## -------------------------------------------------------------------------

pdf("correlation_male_queen_reads.pdf")
getCorrelation(subset_methBase,plot=TRUE)
dev.off()

pdf("sample_cluster_male_queen_reads.pdf")
clusterSamples(subset_methBase, dist="correlation", method="ward", plot=TRUE)
dev.off()

pdf("pca_plot_male_queen_reads.pdf")
PCASamples(subset_methBase)
dev.off()

## -------------------------------------------------------------------------

# Diff meth between worker male and queen read sets
covariates <- data.frame(colony=c("08","19","23","37","08","19","23","37"))
diff_meth <- calculateDiffMeth(subset_methBase, covariates=covariates, mc.cores = 3)
write.csv(diff_meth, file="all_results_diff_meth.csv")

diff_meth_5 <- getMethylDiff(diff_meth, difference=5, qvalue=0.05)
write.csv(diff_meth_5, file="DMRs_min5percentDiff_qval0.05.csv")

# Try no covariates (same result just more significant q-vals)
diff_meth <- calculateDiffMeth(subset_methBase, mc.cores = 3)
write.csv(diff_meth, file="all_results_diff_meth_no_covariates.csv")

diff_meth_5 <- getMethylDiff(diff_meth, difference=5, qvalue=0.05)
write.csv(diff_meth_5, file="DMRs_min5percentDiff_qval0.05_no_covariates.csv")

