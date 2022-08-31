#---------------------------------------------------
# Background gene sets for levels of meth GO analysis
#---------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/GO_analysis")
library(readr)

Bumble_bee_ensemble_GO_terms <- read_delim("Bumble_bee_ensemble_GO_terms.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           col_names = FALSE, trim_ws = TRUE)
colnames(Bumble_bee_ensemble_GO_terms) <- c("gene_id", "GOIds")

#---------------------------------------------------
# Look at genes with high methylation compared to methylated genes per sex/caste
#---------------------------------------------------

# methylated genes background sets
male_all_methylated_genes <- read_csv("male_all_methylated_genes.txt")
queen_all_methylated_genes <- read_csv("queen_all_methylated_genes.txt")
worker_all_methylated_genes <- read_csv("worker_all_methylated_genes.txt")
# Male 7406
# Queen 7944
# Worker 7382

male_all_methylated_genes_GO <- merge(male_all_methylated_genes, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(male_all_methylated_genes_GO$gene_id)) #6319/7406

queen_all_methylated_genes_GO <- merge(queen_all_methylated_genes, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(queen_all_methylated_genes_GO$gene_id)) # 6669/7944

worker_all_methylated_genes_GO <- merge(worker_all_methylated_genes, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(worker_all_methylated_genes_GO$gene_id)) #6318/7382

write.table(male_all_methylated_genes_GO, file="GO_male_all_methylated_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(queen_all_methylated_genes_GO, file="GO_queen_all_methylated_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(worker_all_methylated_genes_GO, file="GO_worker_all_methylated_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)

#---------------------------------------------------
# Differential analysis
#---------------------------------------------------

# Full diff meth lists for background
MvW <- read_delim("MvW_all_diff_meth_genes.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
MvQ <- read_delim("MvQ_all_diff_meth_genes.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
WvQ <- read_delim("WvQ_all_diff_meth_genes.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
# MvW 161
# MvQ 161
# WvQ 59

MvW_GO <- merge(MvW, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(MvW_GO$gene_id)) #151/161

MvQ_GO <- merge(MvQ, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(MvQ_GO$gene_id)) #148/161

WvQ_GO <- merge(WvQ, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(WvQ_GO$gene_id)) #55/59

write.table(MvW_GO, file="GO_diff_meth_MvW.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(MvQ_GO, file="GO_diff_meth_MvQ.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(WvQ_GO, file="GO_diff_meth_WvQ.txt", sep="\t",quote = F, col.names = T, row.names = F)

# Also make a methylated gene list for each comparison for the full diff meth GO
head(male_all_methylated_genes_GO)

MvW_all_background <- rbind(male_all_methylated_genes_GO, worker_all_methylated_genes_GO)
MvW_all_background <- MvW_all_background[!duplicated(MvW_all_background),]

MvQ_all_background <- rbind(male_all_methylated_genes_GO, queen_all_methylated_genes_GO)
MvQ_all_background <- MvQ_all_background[!duplicated(MvQ_all_background),]

WvQ_all_background <- rbind(worker_all_methylated_genes_GO, queen_all_methylated_genes_GO)
WvQ_all_background <- WvQ_all_background[!duplicated(WvQ_all_background),]

write.table(MvW_all_background, file="MvW_all_background.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(MvQ_all_background, file="MvQ_all_background.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(WvQ_all_background, file="WvQ_all_background.txt", sep="\t",quote = F, col.names = T, row.names = F)


# Also want to do hypermeth in Q+W unique to male
all_diff_meth_geneIDs_with_category <- read_delim("all_diff_meth_geneIDs_with_category.txt", 
                                                  delim = "\t", escape_double = FALSE, 
                                                  trim_ws = TRUE)
head(all_diff_meth_geneIDs_with_category)

all_comparisons <- as.data.frame(all_diff_meth_geneIDs_with_category$gene_id)
colnames(all_comparisons) <- "gene_id"
all_comparisons_GO <- merge(all_comparisons, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(all_comparisons_GO$gene_id)) # 188/205

write.table(all_comparisons_GO, file="all_comparisons_background.txt", sep="\t",quote = F, col.names = T, row.names = F)

# Hyper in worker and queen compared to male n = 28
hyper_WQvM <- as.data.frame(all_diff_meth_geneIDs_with_category$gene_id[all_diff_meth_geneIDs_with_category$`Hypermethylated in workers compared to males` ==1 &
                                                    all_diff_meth_geneIDs_with_category$`Hypermethylated in queens compared to males` ==1])
colnames(hyper_WQvM) <- "gene_id"
write.table(hyper_WQvM, file="GO_diff_meth_hyper_WQvM.txt", sep="\t",quote = F, col.names = T, row.names = F)

# Hyper in male compared to worker and queen n = 23
hyper_MvWQ <- as.data.frame(all_diff_meth_geneIDs_with_category$gene_id[all_diff_meth_geneIDs_with_category$`Hypermethylated in males compared to workers` ==1 &
                                                                          all_diff_meth_geneIDs_with_category$`Hypermethylated in males compared to queens` ==1])
colnames(hyper_MvWQ) <- "gene_id"
write.table(hyper_MvWQ, file="GO_diff_meth_hyper_MvWQ.txt", sep="\t",quote = F, col.names = T, row.names = F)
