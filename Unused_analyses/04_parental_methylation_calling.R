#------------------------------------------------
# Calling parent-specific methylation
#------------------------------------------------

# The input alignments here are the worker alignments to either the paretnal male or queen 
# i.e. m08 and q08 are the same worker sample but aligned to either the father or mother of those pooled individuals

library(methylKit)
library(readr)

# Process bam files to make easier to use text files for trialling multiple runs of this script
file.list <- list("w08_male_reads_sorted.bam","w19_male_reads_sorted.bam",
                  "w23_male_reads_sorted.bam","w37_male_reads_sorted.bam",
                  "w08_queen_reads_sorted.bam","w19_queen_reads_sorted.bam", 
                  "w23_queen_reads_sorted.bam","w37_queen_reads_sorted.bam")

raw_data <- processBismarkAln(file.list,
                              sample.id = list("w08_male","w19_male","w23_male","w37_male",
                                               "w08_queen","w19_queen", "w23_queen", "w37_queen"),
                              treatment = c(1,1,1,1,0,0,0,0),
                              assembly="bter_1.0", 
                              read.context="CpG",
                              save.context = "CpG",
                              save.folder=getwd(),
                              save.db = TRUE)


# Read in previously created text files
file.list <- list("w08_male_CpG.txt", "w19_male_CpG.txt", "w23_male_CpG.txt", "w37_male_CpG.txt", 
                  "w08_queen_CpG.txt", "w19_queen_CpG.txt", "w23_queen_CpG.txt", "w37_queen_CpG.txt")

raw_data <- methRead(file.list,
                     sample.id = list("w08_male","w19_male","w23_male","w37_male",
                                      "w08_queen","w19_queen", "w23_queen", "w37_queen"),
                     treatment = c(1,1,1,1,0,0,0,0),
                     assembly="bter_1.0", 
                     context="CpG")


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
                                      "w23_queen", "w37_queen"),
                     treatment = c(1,1,0,0),
                     assembly="bter_1.0", 
                     context="CpG")


# Filter by coverage NOTE: check how much data we lose here, remember the reads are split between two genomes
filtered_data <- filterByCoverage(raw_data,lo.count=8,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)

# Select only CpGs found in all alignments
meth_all_data <- unite(filtered_data, destrand=TRUE) 
nrow(meth_all_data) 
# 24 for all data (N genome)
# cross 1:1218 (N genome)
# cross 2:767 (N genome)

# all data (alt genome): 4
# cross 1 (alt genome): 864
# cross 2 (alt genome): 503

# new genome, new SNPs, N-masked:
# all:
# cross 1:
# cross 2:

## -------------------------------------------------------------------------

# Filter sites using a binomial test so only keep CpGs which are methylated in at least one sample
df_meth_all <- getData(meth_all_data)

a <- df_meth_all[,5:6]
b <- df_meth_all[,8:9]
c <- df_meth_all[,11:12]
d <- df_meth_all[,14:15]
e <- df_meth_all[,17:18]
f <- df_meth_all[,20:21]
g <- df_meth_all[,23:24]
h <- df_meth_all[,26:27]

# NOTE: p shouold be the average non-conversion rate (proportion of methylated Cs compared to non-meth Cs)
# So if 1000 methylated Cs compared to 200,000 T's then 1000/200,000 = 0.005
# for a paper: 'the success probability is the non-conversion rate'
bt <- function(a, b, p = 0.005) {binom.test(a, b, 0.005, alternative="greater")$p.value}
allrows <- data.frame(CT=numeric(),Ccount=numeric(),pVal=numeric(),FDR=numeric(),row=numeric())

for (df in list(a,b,c,d,e,f,g,h)) {
#for (df in list(a,b,c,d)) {
  colnames(df) <- c("CT", "Ccount")
  df$pVal <- mapply(bt, df$Ccount, df$CT)
  df$FDR <- p.adjust(df$pVal, method = "BH", n = length(df$pVal))
  df$row <- as.numeric(rownames(df))
  dfmeth <- subset(df, df$FDR < 0.05)
  allrows <- rbind(allrows, dfmeth)
}
meth_positions <- as.vector(as.numeric(unique(allrows$row))) 
length(meth_positions) 
# cross 2 (N genome): 0
# cross 1 (N genome): 0

# cross 2 (alt genome): 3
# cross 1 (alt genome): 1

# new genome, new SNPs, N-masked:
# all:
# cross 1:
# cross 2:

subset_methBase <- methylKit::select(meth_all_data, meth_positions)

# eyeballed the 3 and 1 and there will be no significant difference 
# for either lineage or parent of origin




pdf("correlation_male_queen.pdf")
getCorrelation(subset_methBase,plot=TRUE)
dev.off()

pdf("sample_cluster_male_queen.pdf")
clusterSamples(subset_methBase, dist="correlation", method="ward", plot=TRUE)
dev.off()

pdf("pca_plot_male_queen.pdf")
PCASamples(subset_methBase)
dev.off()


# Diff meth between worker alignments to parental genomes
covariates <- data.frame(colony=c("08","19","23","37","08","19","23","37"))
diff_meth <- calculateDiffMeth(subset_methBase, covariates=covariates, mc.cores = 1)
write.csv(diff_meth, file="all_results_diff_meth.csv")

diff_meth_5 <- getMethylDiff(diff_meth, difference=5, qvalue=0.05)
write.csv(diff_meth_5, file="DMRs_min5percentDiff_qval0.05_MSCfilter.csv")

# Try no covariates (same result just more significant q-vals)
diff_meth <- calculateDiffMeth(subset_methBase, mc.cores = 3)
write.csv(diff_meth, file="all_results_diff_meth_no_covariates.csv")

diff_meth_5 <- getMethylDiff(diff_meth, difference=5, qvalue=0.05)
write.csv(diff_meth_5, file="DMRs_min5percentDiff_qval0.05_MSCfilter_no_covariates.csv")











#write.table(data_for_model, file= "methylated_positions_with_counts_from_methylkit.txt",
#            sep="\t", col.names=T, row.names=F, quote=F)


# Merge rows to get counts per male alignment and female alignments



# Add in additional columns for the model
#data_for_model$familyID="01-02" 
#data_for_model$familyID[data_for_model$sample="m23"|data_for_model$sample="q37"
#                        |data_for_model$sample="m37"|data_for_model$sample="q23"]="03-04"

#data_for_model$direction_cross="initial"
#data_for_model$direction_cross[data_for_model$sample="m19"|data_for_model$sample="m37"
 #                              |data_for_model$sample="q19"|data_for_model$sample="q37"]="reciprocal"

#data_for_model$BS_percentage_error_rate="0.5"

# Write out dataframe for use in the parent-of-origin methylation model
#write.table(data_for_model, file= "methylated_positions_with_counts_from_methylkit.txt",
#            sep="\t", col.names=T, row.names=F, quote=F)












