# --------------------------------------------------------
# Find out the number of possible syn and non-syn sites per gene
# --------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022")

library(readr)
library(doBy)
library(tidyr)

# Make the gene to transcript matching up file
#grep "mRNA" ref_Bter_1.0_top_level.gff3 > mrna
#cut -f9 mrna > mrna2
#sed -i 's/.*gene=//g' mrna2 
#sed -i 's/;.*transcript_id=/ /g' mrna2
#sed -i 's/;.*protein_id=/ /g' mrna2
#sed 's/;.*//g' mrna2 > gene_transcript.txt

degenerate_sites <- read_delim("Bter_0x4x.txt", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)
tail(degenerate_sites, n=100)

# cut it down
degenerate_sites <- degenerate_sites[,c(5,10)]

# Add synonymous count
degenerate_sites$syn <- 1
degenerate_sites$syn[degenerate_sites$degeneracy==3] <- 0.66
degenerate_sites$syn[degenerate_sites$degeneracy==2] <- 0.5
degenerate_sites$syn[degenerate_sites$degeneracy==0] <- 0

# Add non-syn count
degenerate_sites$non_syn <- 0
degenerate_sites$non_syn[degenerate_sites$degeneracy==3] <- 0.33
degenerate_sites$non_syn[degenerate_sites$degeneracy==2] <- 0.5
degenerate_sites$non_syn[degenerate_sites$degeneracy==0] <- 1

# summary by transcript
degenerate_sites <- degenerate_sites[,-2]
summary_data <- summaryBy(syn + non_syn ~ transcript, data=degenerate_sites, FUN=sum) 
head(summary_data)

# Add in gene names to match the transcript names
gene_transcript <- read_table2("./degenerate_sites/gene_transcript.txt", 
                               col_names = FALSE)
colnames(gene_transcript)<-c("gene_id","transcript")

summary_data <- merge(summary_data, gene_transcript, by="transcript")
summary_data <- summary_data[complete.cases(summary_data),]

# summary by gene id
summary_data <- summary_data[,-1]
summary_data <- summaryBy(syn.sum + non_syn.sum ~ gene_id, data=summary_data, FUN=sum) 
colnames(summary_data)<-c("gene_id","syn_all_gene","non_syn_all_gene")

head(summary_data)
write.table(summary_data, file="avaliable_syn_nonsyn.txt", sep="\t", quote=F, col.names = T, row.names = F)
