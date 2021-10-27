#------------------------------------------------
# Calling lineage-specific methylation
#------------------------------------------------
library(methylKit)
library(readr)
setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/Honeybee_files/methylkit_inputs")

# --- for each cross

sample.list <- list("875_1_drone.CpG_report.merged_CpG_evidence.cov",
                    "875_1_queen.CpG_report.merged_CpG_evidence.cov",
                    "875_2_drone.CpG_report.merged_CpG_evidence.cov",
                    "875_2_queen.CpG_report.merged_CpG_evidence.cov",
                    "875_E_drone.CpG_report.merged_CpG_evidence.cov",
                    "875_E_queen.CpG_report.merged_CpG_evidence.cov",
                    "875_F_drone.CpG_report.merged_CpG_evidence.cov",
                    "875_F_queen.CpG_report.merged_CpG_evidence.cov",
                    "888_1_drone.CpG_report.merged_CpG_evidence.cov",
                    "888_1_queen.CpG_report.merged_CpG_evidence.cov",
                    "888_2_drone.CpG_report.merged_CpG_evidence.cov",
                    "888_2_queen.CpG_report.merged_CpG_evidence.cov",
                    "888_E_drone.CpG_report.merged_CpG_evidence.cov",
                    "888_E_queen.CpG_report.merged_CpG_evidence.cov",
                    "888_F_drone.CpG_report.merged_CpG_evidence.cov",
                    "888_F_queen.CpG_report.merged_CpG_evidence.cov")

CPGRaw <- methRead(sample.list, 
                   sample.id = list("875_1_drone",
                                    "875_1_queen",
                                    "875_2_drone",
                                    "875_2_queen",
                                    "875_E_drone",
                                    "875_E_queen",
                                    "875_F_drone",
                                    "875_F_queen",
                                    "888_1_drone",
                                    "888_1_queen",
                                    "888_2_drone",
                                    "888_2_queen",
                                    "888_E_drone",
                                    "888_E_queen",
                                    "888_F_drone",
                                    "888_F_queen"),
                   assembly="apis3",
                   treatment=c(0,0,0,0,0,0,0,0,
                               1,1,1,1,1,1,1,1),
                   context="CpG",
                   dbtype = NA,
                   pipeline = "bismarkCoverage",
                   header = T, sep = "\t", mincov=1,
                   dbdir= path)


sample.list <- list("882_0_drone.CpG_report.merged_CpG_evidence.cov",
                    "882_0_queen.CpG_report.merged_CpG_evidence.cov",
                    "882_2_drone.CpG_report.merged_CpG_evidence.cov",
                    "882_2_queen.CpG_report.merged_CpG_evidence.cov",
                    "882_3_drone.CpG_report.merged_CpG_evidence.cov",
                    "882_3_queen.CpG_report.merged_CpG_evidence.cov",
                    "882_5_drone.CpG_report.merged_CpG_evidence.cov",
                    "882_5_queen.CpG_report.merged_CpG_evidence.cov",
                    "894_1_drone.CpG_report.merged_CpG_evidence.cov",
                    "894_1_queen.CpG_report.merged_CpG_evidence.cov",
                    "894_2_drone.CpG_report.merged_CpG_evidence.cov",
                    "894_2_queen.CpG_report.merged_CpG_evidence.cov",
                    "894_5_drone.CpG_report.merged_CpG_evidence.cov",
                    "894_5_queen.CpG_report.merged_CpG_evidence.cov",
                    "894_6_drone.CpG_report.merged_CpG_evidence.cov",
                    "894_6_queen.CpG_report.merged_CpG_evidence.cov")

CPGRaw <- methRead(sample.list, 
                   sample.id = list("882_0_drone",
                                    "882_0_queen",
                                    "882_2_drone",
                                    "882_2_queen",
                                    "882_3_drone",
                                    "882_3_queen",
                                    "882_5_drone",
                                    "882_5_queen",
                                    "894_1_drone",
                                    "894_1_queen",
                                    "894_2_drone",
                                    "894_2_queen",
                                    "894_5_drone",
                                    "894_5_queen",
                                    "894_6_drone",
                                    "894_6_queen"),
                   assembly="apis3",
                   treatment=c(0,0,0,0,0,0,0,0,
                               1,1,1,1,1,1,1,1),
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
# cross 1: 6383 (with 4/8 samples = 567336)
# cross 2: 4464 (with 4/8 samples = 360535)

## -------------------------------------------------------------------------

# Filter sites using a binomial test so only keep CpGs which are methylated in at least one sample
df_meth_all <- getData(meth_all_data)

# Change all NAs to something so the test runs, this is ok as long as we don't change them to 
# something which calls as methylated, as the function of the binomial here is to select
# sites which are methylated in at least one sample, so doing this won't affect the sites we select
# but it will allow the binomial to run and maintain the correct number of rows per sample, this is
# super important as we use the row numbers to later select the correct CpGs from the methylkit object

# This also works because we only run the binomial test to get the row numbers we want to keep
# we then subset the methylkit object using these row numbers so when we run the differential methylation
# test we are still running it on the original data which shows NA when appropriate, not the new values

# Ugly way of changing every NA coverage value to 100
df_meth_all[,c(5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50)][is.na
            (df_meth_all[,c(5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50)])] <- 100
# Ugly way of changing every numCs value to 1
df_meth_all[,c(6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51)][is.na
            (df_meth_all[,c(6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51)])] <- 1
# Using these values means the NA site is not called as methylated
# This is a lazy way of doing it and not sensible if you have 100s of samples

a <- df_meth_all[,5:6]
b <- df_meth_all[,8:9]
c <- df_meth_all[,11:12]
d <- df_meth_all[,14:15]
e <- df_meth_all[,17:18]
f <- df_meth_all[,20:21]
g <- df_meth_all[,23:24]
h <- df_meth_all[,26:27]

k <- df_meth_all[,29:30]
l <- df_meth_all[,32:33]
m <- df_meth_all[,35:36]
n <- df_meth_all[,38:39]
o <- df_meth_all[,41:42]
p <- df_meth_all[,44:45]
q <- df_meth_all[,47:48]
r <- df_meth_all[,50:51]


# NOTE: p shouold be the average non-conversion rate (proportion of methylated Cs compared to non-meth Cs)
# So if 1000 methylated Cs compared to 200,000 T's then 1000/200,000 = 0.005
# for a paper: 'the success probability is the non-conversion rate'
bt <- function(a, b, p = 0.005) {binom.test(a, b, 0.005, alternative="greater")$p.value}
allrows <- data.frame(CT=numeric(),Ccount=numeric(),pVal=numeric(),FDR=numeric(),row=numeric())

for (df in list(a,b,c,d,e,f,g,h,k,l,m,n,o,p,q,r)) {
  colnames(df) <- c("CT", "Ccount")
  df$pVal <- mapply(bt, df$Ccount, df$CT)
  df$FDR <- p.adjust(df$pVal, method = "BH", n = length(df$pVal))
  df$row <- as.numeric(rownames(df))
  dfmeth <- subset(df, df$FDR < 0.05)
  allrows <- rbind(allrows, dfmeth)
}
meth_positions <- as.vector(as.numeric(unique(allrows$row))) 
length(meth_positions) 
# cross 1: 26
# cross 2: 22


subset_methBase <- methylKit::select(meth_all_data, meth_positions)
head(subset_methBase)

# Diff meth between worker alignments to parental genomes
covariates <- data.frame(colony=c("male","female","male","female","male","female","male","female",
                                  "male","female","male","female","male","female","male","female"))

diff_meth <- calculateDiffMeth(subset_methBase, covariates=covariates, mc.cores = 1)
write.csv(diff_meth, file="all_results_diff_meth.csv")

diff_meth_5 <- getMethylDiff(diff_meth, difference=5, qvalue=0.05)
write.csv(diff_meth_5, file="DMRs_min5percentDiff_qval0.05_MSCfilter.csv")

# cross 1: 16/26
# cross 2: 6/22