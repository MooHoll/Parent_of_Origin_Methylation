# Take a look at the diff meth sites identified on ALICE for leuven-meth project

library(readr)

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Methylation/Diff_Meth")

# Worker v males DMRs
w_m <- read_csv("workers_vs_males/DMRs_min10percentDiff_qval0.05_MSCfilter_worker_male.csv")
w_m_sig <- subset(w_m, pvalue<0.05) #none

# Worker v queens DMRs
w_q <- read_csv("queen_vs_workers/DMRs_min10percentDiff_qval0.05_MSCfilter_worker_queen.csv")
w_q_sig <- subset(w_q, pvalue<0.05) #none

# male v queens DMRs
q_m <- read_csv("queens_vs_males/DMRs_min10percentDiff_qval0.05_MSCfilter_worker_queen.csv")
q_m_sig <- subset(q_m, qvalue<0.05) #703 cool! but what about ploidy???

### ON ALICE

library(sqldf)
library(readr)

genome_annotation<-read.csv.sql("just_genes.csv",
                                sql ="select * from file", sep=",",header = T)
diff_meth_sites <- read.csv("queens_vs_males/DMRs_min10percentDiff_qval0.05_MSCfilter_worker_queen.csv")
colnames(diff_meth_sites)[8]<-"meth_diff"

output <- sqldf("SELECT diff_meth_sites.chr,
      diff_meth_sites.start,
                diff_meth_sites.end,
                diff_meth_sites.strand,
                diff_meth_sites.pvalue,
                diff_meth_sites.qvalue,
                diff_meth_sites.meth_diff,
                ga.chr,
                ga.start,
                ga.end,
                ga.geneID
                FROM diff_meth_sites AS diff_meth_sites
                LEFT JOIN genome_annotation AS ga 
                ON diff_meth_sites.chr = ga.chr
                AND (diff_meth_sites.start >= ga.start AND diff_meth_sites.start <= ga.end)") 


output_subset_no_geneIDs <-subset(output, !output$geneID=="NA") #660 out of 703 had gene IDs
output_dedup <- subset(output_subset_no_geneIDs, !duplicated(geneID)) # of those 412 unique genes

write.table(output_subset_no_geneIDs, file="DMGs_MSCfilter_with_geneID_Male_Queen.csv", sep="\t", row.names=F, quote=F)

write.csv(output_subset_no_geneIDs$geneID, file="list_DMGs_MSCfilter_Male_Queen.csv")



