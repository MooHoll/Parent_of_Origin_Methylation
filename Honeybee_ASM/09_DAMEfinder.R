# ------------------------------------------------------------------------
# Identification of differential ASM
# ------------------------------------------------------------------------

# Load packages
library(DAMEfinder)
library(SummarizedExperiment)
library(dplyr)

# Set working directory
setwd("~/Dropbox/Leicester_postdoc/Students/Developmental_ASM_2021/Data/Bumblebee")

# Define which files we want
tuple_files <- c("A1_sorted_methtuple.CG.2.tsv.gz","A2_sorted_methtuple.CG.2.tsv.gz",
                         "A3_sorted_methtuple.CG.2.tsv.gz","A4_sorted_methtuple.CG.2.tsv.gz",
                         "B1_sorted_methtuple.CG.2.tsv.gz","B2_sorted_methtuple.CG.2.tsv.gz",
                         "B3_sorted_methtuple.CG.2.tsv.gz","B4_sorted_methtuple.CG.2.tsv.gz")
# Give them meaningful names
sample_names <- c("repro_brain_rep1","repro_brain_rep2","repro_brain_rep3","repro_brain_rep4",
                  "repro_ovary_rep1","repro_ovary_rep2","repro_ovary_rep3","repro_ovary_rep4")

# Read in the files
tuple_list <- read_tuples(tuple_files, sample_names)

# Before we calculate the ASM scores we want to filter to ensure we have at least 10 coverage
for (i in seq_along(tuple_list)){
  tuple_list[[i]] <- tuple_list[[i]][tuple_list[[i]]$cov>9,]
}

# Calculate the ASM scores
ASMscore <- calc_asm(tuple_list)

# Filter out none-calls
filter_none_calls <- rowSums(!is.na(assay(ASMscore, "asm"))) == 8
ASMscore <- ASMscore[filter_none_calls,]

# Make design matrix (or specify a contrast)
grp <- factor(c(rep("repro_brain",4),rep("repro_ovary",4)), levels = c("repro_brain", "repro_ovary"))
mod <- model.matrix(~grp)

# Make an MDS plot of the samples
methyl_MDS_plot(ASMscore, group = grp)

# Filter data to keep only ASM regions which score > 0.8 in at least one sample
# NOTE: 0.8 is the score deemed to be representative of true ASM according to the
# DAMEfinder publication
filter_0.8 <- as.data.frame(SummarizedExperiment::assays(ASMscore)$asm)
filter_0.8 <-  filter_0.8 %>% filter_all(any_vars(.>0.8))
filter_0.8$positions <- row.names(filter_0.8)

ASMscore_subset <- ASMscore[filter_0.8$positions,]

# Run in default mode and specify the max gap between CpGs
dames <- find_dames(ASMscore_subset, mod, maxGap = 100)

# Filter to keep only significant dames
sig_dames <- dames[dames$FDR < 0.05,]
write.table(sig_dames, file="repro_brain_vs_repro_ovary_DAMES.txt", sep="\t", 
            quote = F, col.names = T, row.names = F)

# Have a look at one DAME as an example, first set names for plotting purposes
colData(ASMscore)$group <- grp
colData(ASMscore)$samples <- colnames(ASMscore)

# Choose region to look at based on the chromosome and genomic coordinates
dame <- GRanges("NW_003570989",IRanges(374,464))
dame_track(dame = dame,ASM = ASMscore)
dame_track_mean(dame = dame,ASM = ASMscore)


