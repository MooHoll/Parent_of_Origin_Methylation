#------------------------------------------------
# Differential methylation between castes
#------------------------------------------------
setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/coverage_files")
library(methylKit)
library(readr)


sample.list <- list("m08_merged_CpG_evidence.cov" ,"m19_merged_CpG_evidence.cov",
                    "m23_merged_CpG_evidence.cov","m37_merged_CpG_evidence.cov",
                    "w08_merged_CpG_evidence.cov" ,"w19_merged_CpG_evidence.cov",
                    "w23_merged_CpG_evidence.cov","w37_merged_CpG_evidence.cov")

CPGRaw <- methRead(sample.list, 
                   sample.id = list("m08", "m19","m23","m37",
                                    "w08", "w19","w23","w37"),
                   assembly="bter_1.0",
                   treatment=c(0,0,0,0,1,1,1,1),
                   context="CpG",
                   dbtype = NA,
                   pipeline = "bismarkCoverage",
                   header = T, sep = "\t", mincov=1,
                   dbdir= path)

# ---

sample.list <- list("m08_merged_CpG_evidence.cov" ,"m19_merged_CpG_evidence.cov",
                    "m23_merged_CpG_evidence.cov","m37_merged_CpG_evidence.cov",
                    "q08_merged_CpG_evidence.cov" ,"q19_merged_CpG_evidence.cov",
                    "q23_merged_CpG_evidence.cov","q37_merged_CpG_evidence.cov")

CPGRaw <- methRead(sample.list, 
                   sample.id = list("m08", "m19","m23","m37",
                                    "q08", "q19","q23","q37"),
                   assembly="bter_1.0",
                   treatment=c(0,0,0,0,1,1,1,1),
                   context="CpG",
                   dbtype = NA,
                   pipeline = "bismarkCoverage",
                   header = T, sep = "\t", mincov=1,
                   dbdir= path)

# ---

sample.list <- list("w08_merged_CpG_evidence.cov" ,"w19_merged_CpG_evidence.cov",
                    "w23_merged_CpG_evidence.cov","w37_merged_CpG_evidence.cov",
                    "q08_merged_CpG_evidence.cov" ,"q19_merged_CpG_evidence.cov",
                    "q23_merged_CpG_evidence.cov","q37_merged_CpG_evidence.cov")

CPGRaw <- methRead(sample.list, 
                   sample.id = list("w08", "w19","w23","w37",
                                    "q08", "q19","q23","q37"),
                   assembly="bter_1.0",
                   treatment=c(0,0,0,0,1,1,1,1),
                   context="CpG",
                   dbtype = NA,
                   pipeline = "bismarkCoverage",
                   header = T, sep = "\t", mincov=1,
                   dbdir= path)

# Filter by coverage
filtered_data <- filterByCoverage(CPGRaw,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)

normalized <- normalizeCoverage(filtered_data)

# Select only CpGs found in all alignments
meth_all_data <- unite(normalized, destrand=F) 
nrow(meth_all_data) 
# male vs worker: 2,448,114
# male vs queen: 2,253,693
# worker vs queen: 2,185,909

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
  colnames(df) <- c("CT", "Ccount")
  df$pVal <- mapply(bt, df$Ccount, df$CT)
  df$FDR <- p.adjust(df$pVal, method = "BH", n = length(df$pVal))
  df$row <- as.numeric(rownames(df))
  dfmeth <- subset(df, df$FDR < 0.05)
  allrows <- rbind(allrows, dfmeth)
}

meth_positions <- as.vector(as.numeric(unique(allrows$row))) 
length(meth_positions) 
# male vs worker: 9298
# male vs queen: 7373
# worker vs queen: 9519

subset_methBase <- methylKit::select(meth_all_data, meth_positions)
head(subset_methBase)

# Diff meth between worker alignments to parental genomes
covariates1 <- data.frame(subspecies=c("audax","dalmatinus", "audax", "dalmatinus", 
                                      "mixed", "mixed", "mixed", "mixed"))

#covariates <- data.frame(subspecies=c("audax","dalmatinus", "audax", "dalmatinus", 
#                                      "audax","dalmatinus", "audax", "dalmatinus"))

#covariates <- data.frame(subspecies=c("audax","dalmatinus", "audax", "dalmatinus", 
#                                      "mixed", "mixed", "mixed", "mixed"))

diff_meth <- calculateDiffMeth(subset_methBase, mc.cores = 1, covariates = covariates1, method='qvalue')

diff_meth_5_covar <- getMethylDiff(diff_meth, difference=10, qvalue=0.01)
nrow(diff_meth_5_covar)

#write.csv(diff_meth_5, file="male_vs_worker_diff_meth_CpGs.csv")
#write.csv(diff_meth_5, file="male_vs_queen_diff_meth_CpGs.csv")
write.csv(diff_meth_5, file="worker_vs_queen_diff_meth_CpGs.csv")

# male vs worker (+ve = worker hyper): 1232
# male vs queen (+ve = queen hyper):  1034
# worker vs queen (+ve = queen hyper): 358