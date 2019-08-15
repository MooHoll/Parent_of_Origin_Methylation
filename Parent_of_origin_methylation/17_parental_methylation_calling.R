#------------------------------------------------
# Calling parent-specific methylation
#------------------------------------------------

# The input alignments here are the worker alignments to either the paretnal male or queen 
# i.e. m08 and q08 are the same worker sample but aligned to either the father or mother of those pooled individuals

library(methylKit)
library(readr)

# Input sam filess must be sorted with fancy command, see script 16.
file.list <- list("m08_no_mismatches_deduplicated_sorted.sam","m19_no_mismatches_deduplicated_sorted.sam",
                  "m23_no_mismatches_deduplicated_sorted.sam","m37_no_mismatches_deduplicated_sorted.sam",
                  "q08_no_mismatches_deduplicated_sorted.sam","q19_no_mismatches_deduplicated_sorted.sam", 
                  "q23_no_mismatches_deduplicated_sorted.sam","q37_no_mismatches_deduplicated_sorted.sam")

raw_data <- processBismarkAln(file.list,
                              sample.id = list("m08","m19","m23","m37","q08","q19", "q23", "q37"),
                              treatment = c(1,1,1,1,0,0,0,0),
                              assembly="bter_1.0", 
                              read.context="CpG",
                              save.context = "CpG",
                              save.folder=getwd(),
                              save.db = TRUE)

# Read in previously created text files
file.list <- list("m08_CpG.txt", "m19_CpG.txt", "m23_CpG.txt", "m37_CpG.txt", 
                  "q08_CpG.txt", "q19_CpG.txt", "q23_CpG.txt", "q37_CpG.txt")

raw_data <- methRead(file.list,
                     sample.id = list("m08","m19","m23","m37","q08","q19", "q23", "q37"),
                     treatment = c(1,1,1,1,0,0,0,0),
                     assembly="bter_1.0", 
                     context="CpG")

# Filter by coverage NOTE: check how much data we lose here, remember the reads are split between two genomes
filtered_data <- filterByCoverage(raw_data,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)

# Select only CpGs found in all alignments
meth_all_data <- unite(filtered_data, destrand=TRUE)
df_meth_all <- getData(meth_all_data)
df_meth_all$rownums <- row.names(df_meth_all)

source("/scratch/monoallelic/hm257/MSC_scripts/MSC.R")
source("/scratch/monoallelic/hm257/MSC_scripts/rateestimate.R")

# Determine positions shown to be methylated per sample
m08 <- df_meth_all[,5:6]
m19 <- df_meth_all[,8:9]
m23 <- df_meth_all[,11:12]
m37 <- df_meth_all[,14:15]
q08 <- df_meth_all[,17:18]
q19 <- df_meth_all[,20:21]
q23 <- df_meth_all[,23:24]
q37 <- df_meth_all[,26:27]

samples <- list(m08, m19, m23, m37, q08, q19, q23, q37)

for(i in seq_along(samples)){
  colnames(samples[[i]]) <- c("CT","Ccount")
  MSCount <- MSC(samples[[i]], 1e-08)
  MSCresult <- MSCount$MSC
  myfile <- file.path("./", paste0(i,"_","msc_meth_calls.csv"))
  write.csv(MSCresult, file=myfile)
  pi <- MSCount$pi
  MSCrate <- rateestimate(MSCresult,pi)
  rate <- as.data.frame(MSCrate)
  myfile2 <- file.path("./", paste0(i,"_","msc_rate_extimates.csv"))
  write.csv(rate, file=myfile2)
  MSCresult$row_nums <- row.names(MSCresult)
  label <- paste("MSCresult_meth_", i, sep="_")
  assign(label, subset(MSCresult, status == "methylated")) 
}


# Get a dataframe which keeps only positions that show methylation in at least one sample
all_data <- rbind(MSCresult_meth__1, MSCresult_meth__2, MSCresult_meth__3, MSCresult_meth__4,
                  MSCresult_meth__5, MSCresult_meth__6, MSCresult_meth__7, MSCresult_meth__8)

meth_positions <- as.vector(as.numeric(unique(all_data$row_nums))) 

subset_methBase <- select (meth_all_data, meth_positions)

pdf("correlation_male_queen.pdf")
getCorrelation(subset_methBase,plot=TRUE)
dev.off()

pdf("sample_cluster_male_queen.pdf")
clusterSamples(subset_methBase, dist="correlation", method="ward", plot=TRUE)
dev.off()

pdf("pca_plot_male_queen.pdf")
PCASamples(subset_methBase)
dev.off()

# Make output files that have the methylation level of each CpG that was tested
methBase_ob <- getData(subset_methBase)
write.table(methBase_ob, file="objectmethbase.txt", quote=F, 
            row.names = F, sep = '\t')

objectmethbase <- read_delim("objectmethbase.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)

objectmethbase$chrBase <- paste(objectmethbase$chr, ".", objectmethbase$start, sep="")
objectmethbase <- objectmethbase[,-3]
objectmethbase$strand <- gsub(pattern = "\\+", replace = "F", objectmethbase$strand)

m08 <- objectmethbase[,c("chrBase","chr","start","strand","coverage1","numCs1","numTs1")]
m19 <- objectmethbase[,c("chrBase","chr","start","strand","coverage2","numCs2","numTs2")]
m23 <- objectmethbase[,c("chrBase","chr","start","strand","coverage3","numCs3","numTs3")]
m37 <- objectmethbase[,c("chrBase","chr","start","strand","coverage4","numCs4","numTs4")]
q08 <- objectmethbase[,c("chrBase","chr","start","strand","coverage5","numCs5","numTs5")]
q19 <- objectmethbase[,c("chrBase","chr","start","strand","coverage6","numCs6","numTs6")]
q23 <- objectmethbase[,c("chrBase","chr","start","strand","coverage7","numCs7","numTs7")]
q37 <- objectmethbase[,c("chrBase","chr","start","strand","coverage8","numCs8","numTs8")]

all_files <- list(m08, m19, m23, m37, q08, q19, q23, q37)

for(i in seq_along(all_files)){
  colnames(all_files[[i]])[c(3,5,6,7)] <- c("base","coverage","numCs","numTs")
  all_files[[i]]$freqC <- round((all_files[[i]]$numCs/all_files[[i]]$coverage)*100,
                                digits=2)
  all_files[[i]]$freqT <- round((all_files[[i]]$numTs/all_files[[i]]$coverage)*100,
                                digits=2)
  all_files[[i]] <- all_files[[i]][-c(6,7)]
  myfile <- file.path("./", paste0(i,"_","subsetted_final.txt"))
  write.table(all_files[[i]], file=myfile, quote=F, sep="\t", row.names=F)
}


# Diff meth between worker alignments to parental genomes
covariates <- data.frame(colony=c("08","19","23","37","08","19","23","37"))
diff_meth <- calculateDiffMeth(subset_methBase, covariates=covariates, mc.cores = 3)
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












