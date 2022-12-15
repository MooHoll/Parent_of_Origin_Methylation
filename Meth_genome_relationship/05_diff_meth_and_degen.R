#-----------------------------------------------------------------------
# Checking overlap between diff meth genes and degen genes
#-----------------------------------------------------------------------
setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/degenerate_sites")

library(readr)
library(reshape2)
library(dplyr)
library(ggplot2)

#-----------------------------------------------------------------------
# Need to make the degen and meth sample file for males and females on average

# binomial methylation calls
methylation_calls_per_sample <- read_delim("methylation_calls_per_sample.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)
head(methylation_calls_per_sample) # ~2500 CpGs methylated per sample

# genome degenreacy (one at a time as files are huge - like 3Gb each)
degenerate_sites <- read_delim("m08_degenerate_sites.txt.gz", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)
head(degenerate_sites)
degenerate_sites <- degenerate_sites[,-c(2,4)]

# Keep only C/G sites
degenerate_sites <- degenerate_sites[degenerate_sites$codon.base == "C" |
                                       degenerate_sites$codon.base == "G",]

colnames(degenerate_sites)<-c("chr","position","gene_id","codon","base","amino_acid",
                              "codon_position","degeneracy")

# Change this depending on sample
sample <- methylation_calls_per_sample[,c(1,2,3)]

# Merge by chromosome and CpG position
all_data <- merge(sample, degenerate_sites, by=c("chr","position"))
colnames(all_data)[3] <- "sample"

all_data$sample[all_data$sample==0] <- "Unmethylated"
all_data$sample[all_data$sample==1] <- "Methylated"
all_data <- all_data[,c(1,2,3,9)]
colnames(all_data)[3] <- "methylation_status"
all_data <- all_data[!duplicated(all_data),]

# Change depending on sample
all_data$sample <- "m08"
write.table(all_data, file="m08_degen_and_meth.txt", quote=F, col.names =T, row.names = F, sep="\t")

#-----------------------------------------------------------------------

# Read them all back in and make one file
file.list = list.files(("./"),pattern="*_degen_and_meth.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)
samples_all <- as.data.frame(bind_rows(samples))

#-----------------------------------------------------------------------

# Diff meth sites
male_vs_worker_diff_meth_CpGs <- read_csv("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/coverage_files/male_vs_worker_diff_meth_CpGs.csv")
male_vs_worker_diff_meth_CpGs <- male_vs_worker_diff_meth_CpGs[,c(2,3,4,8)]
male_vs_queen_diff_meth_CpGs <- read_csv("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/coverage_files/male_vs_queen_diff_meth_CpGs.csv")
male_vs_queen_diff_meth_CpGs <- male_vs_queen_diff_meth_CpGs[,c(2,3,4,8)]
worker_vs_queen_diff_meth_CpGs <- read_csv("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/coverage_files/worker_vs_queen_diff_meth_CpGs.csv")
worker_vs_queen_diff_meth_CpGs <- worker_vs_queen_diff_meth_CpGs[,c(2,3,4,8)]

# Checked and only the start position relates to the above degen file, no overlap with any end positions
colnames(male_vs_worker_diff_meth_CpGs)[2] <- "position"
male_worker <- merge(samples_all, male_vs_worker_diff_meth_CpGs, by = c("chr", "position"), all.y = T)
colnames(male_vs_queen_diff_meth_CpGs)[2] <- "position"
male_queen <- merge(samples_all, male_vs_queen_diff_meth_CpGs, by = c("chr", "position"), all.y = T)
colnames(worker_vs_queen_diff_meth_CpGs)[2] <- "position"
worker_queen <- merge(samples_all, worker_vs_queen_diff_meth_CpGs, by = c("chr", "position"), all.y = T)

# See how many are even labelled with degeneracy, i.e. in a coding exon
sum(is.na(male_worker$sample)) # 758 / 1232 do not have degeneracy labelled 61%
sum(is.na(male_queen$sample)) # 385 / 1034 = 37%
sum(is.na(worker_queen$sample)) # 187 / 358 = 52%

# Keep only the male queen analysis going forward
head(male_queen)
male_queen <- male_queen[complete.cases(male_queen),]

# Lets have a look if the degen looks particularly different anyway for queens/males
ggplot(male_queen, aes(x=degeneracy, fill=sample))+
  geom_bar(stat="count", position = "dodge") # nope all the same

# So where are the diff meth sites
male_queen$hyper <- "Male"
male_queen$hyper[male_queen$meth.diff > 0] <- "Queen"
male_queen$degeneracy <- as.factor(male_queen$degeneracy)

ggplot(male_queen, aes(x=degeneracy, fill=hyper))+
  geom_bar(stat="count", position = "dodge") 

# This is confounded, of course diff meth are more likely to be in 0x as these are more likely to be meth
# generally, so it's the same thing
head(male_queen)
male_queen <- male_queen[,c(4,7)]
male_queen <- male_queen[!duplicated(male_queen),]

ggplot(male_queen, aes(x=degeneracy))+
  geom_bar(stat="count", position = "dodge")+
  theme_bw()+
  xlab("Codon Degeneracy")+
  ylab("Count")+
  ggtitle("Differentially methylated positions")+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())

# Now need to see if the 0x degen sites occur in highly methylated genes etc
# first question, do any genes even show more of one degeneracy that another,
# or is it spread across genes, i.e. one gene shows 0,2,3,4x

# Ask the basic question first, how many genes have 0x sites



# binomial methylation calls
methylation_calls_per_sample <- read_delim("methylation_calls_per_sample.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)
head(methylation_calls_per_sample) # ~2500 CpGs methylated per sample

# genome degenreacy (one at a time as files are huge - like 3Gb each)
degenerate_sites <- read_delim("m08_degenerate_sites.txt.gz", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)
head(degenerate_sites)
degenerate_sites <- degenerate_sites[,-c(2,4)]

# Keep only C/G sites
degenerate_sites <- degenerate_sites[degenerate_sites$codon.base == "C" |
                                       degenerate_sites$codon.base == "G",]

colnames(degenerate_sites)<-c("chr","position","gene_id","codon","base","amino_acid",
                              "codon_position","degeneracy")

# Change this depending on sample
sample <- methylation_calls_per_sample[,c(1,2,3)]

# Merge by chromosome and CpG position
all_data <- merge(sample, degenerate_sites, by=c("chr","position"))
colnames(all_data)[3] <- "sample"

all_data$sample[all_data$sample==0] <- "Unmethylated"
all_data$sample[all_data$sample==1] <- "Methylated"








