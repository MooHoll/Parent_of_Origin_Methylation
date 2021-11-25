#---------------------------------------------------
# Background gene sets for levels of meth GO analysis
#---------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/GO_analysis")
library(readr)

Bumble_bee_ensemble_GO_terms <- read_delim("Bumble_bee_ensemble_GO_terms.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           col_names = FALSE, trim_ws = TRUE)
colnames(Bumble_bee_ensemble_GO_terms) <- c("gene_id", "GOIds")

#---------------------------------------------------
# Look at genes with high methylation compared to methylated genes per sex/caste
#---------------------------------------------------

# methylated genes background sets
male_all_methylated_genes <- read_csv("meth_level_gene_lists/male_all_methylated_genes.txt")
queen_all_methylated_genes <- read_csv("meth_level_gene_lists/queen_all_methylated_genes.txt")
worker_all_methylated_genes <- read_csv("meth_level_gene_lists/worker_all_methylated_genes.txt")
# Male 7206
# Queen 7754
# Worker 7316

male_all_methylated_genes_GO <- merge(male_all_methylated_genes, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(male_all_methylated_genes_GO$gene_id)) #6139/7206

queen_all_methylated_genes_GO <- merge(queen_all_methylated_genes, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(queen_all_methylated_genes_GO$gene_id)) # 6550/7754

worker_all_methylated_genes_GO <- merge(worker_all_methylated_genes, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(worker_all_methylated_genes_GO$gene_id)) #6252/7316

write.table(male_all_methylated_genes_GO, file="meth_level_gene_lists/GO_male_all_methylated_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(queen_all_methylated_genes_GO, file="meth_level_gene_lists/GO_queen_all_methylated_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(worker_all_methylated_genes_GO, file="meth_level_gene_lists/GO_worker_all_methylated_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)

#---------------------------------------------------
# Differential analysis
#---------------------------------------------------

# Full diff meth lists for background
MvW <- read_delim("diff_meth_gene_lists/MvW_all_diff_meth_genes.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
MvQ <- read_delim("diff_meth_gene_lists/MvQ_all_diff_meth_genes.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
WvQ <- read_delim("diff_meth_gene_lists/WvQ_all_diff_meth_genes.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
# MvW 155
# MvQ 165
# WvQ 37

MvW_GO <- merge(MvW, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(MvW_GO$gene_id)) #146/155

MvQ_GO <- merge(MvQ, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(MvQ_GO$gene_id)) #152/165

WvQ_GO <- merge(WvQ, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(WvQ_GO$gene_id)) #34/37

write.table(MvW_GO, file="diff_meth_gene_lists/GO_diff_meth_MvW.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(MvQ_GO, file="diff_meth_gene_lists/GO_diff_meth_MvQ.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(WvQ_GO, file="diff_meth_gene_lists/GO_diff_meth_WvQ.txt", sep="\t",quote = F, col.names = T, row.names = F)

# Also make a methylated gene list for each comparison for the full diff meth GO
head(male_all_methylated_genes_GO)

MvW_all_background <- rbind(male_all_methylated_genes_GO, worker_all_methylated_genes_GO)
MvW_all_background <- MvW_all_background[!duplicated(MvW_all_background),]

MvQ_all_background <- rbind(male_all_methylated_genes_GO, queen_all_methylated_genes_GO)
MvQ_all_background <- MvQ_all_background[!duplicated(MvQ_all_background),]

WvQ_all_background <- rbind(worker_all_methylated_genes_GO, queen_all_methylated_genes_GO)
WvQ_all_background <- WvQ_all_background[!duplicated(WvQ_all_background),]

write.table(MvW_all_background, file="diff_meth_gene_lists/MvW_all_background.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(MvQ_all_background, file="diff_meth_gene_lists/MvQ_all_background.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(WvQ_all_background, file="diff_meth_gene_lists/WvQ_all_background.txt", sep="\t",quote = F, col.names = T, row.names = F)


# Also want to do hypermeth in Q+W unique to male
all_diff_meth_geneIDs_with_category <- read_delim("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/Differential_methylation/all_diff_meth_geneIDs_with_category.txt", 
                                                  delim = "\t", escape_double = FALSE, 
                                                  trim_ws = TRUE)
head(all_diff_meth_geneIDs_with_category)

all_comparisons <- as.data.frame(all_diff_meth_geneIDs_with_category$gene_id)
colnames(all_comparisons) <- "gene_id"
all_comparisons_GO <- merge(all_comparisons, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(all_comparisons_GO$gene_id)) # 164/178

write.table(all_comparisons_GO, file="diff_meth_gene_lists/all_comparisons_background.txt", sep="\t",quote = F, col.names = T, row.names = F)

# Hyper in worker and queen compared to male n = 31
hyper_WQvM <- as.data.frame(all_diff_meth_geneIDs_with_category$gene_id[all_diff_meth_geneIDs_with_category$`Hypermethylated in workers compared to males` ==1 &
                                                    all_diff_meth_geneIDs_with_category$`Hypermethylated in queens compared to males` ==1])
colnames(hyper_WQvM) <- "gene_id"
write.table(hyper_WQvM, file="diff_meth_gene_lists/GO_diff_meth_hyper_WQvM.txt", sep="\t",quote = F, col.names = T, row.names = F)

# Hyper in male compared to worker and queen n = 18
hyper_MvWQ <- as.data.frame(all_diff_meth_geneIDs_with_category$gene_id[all_diff_meth_geneIDs_with_category$`Hypermethylated in males compared to workers` ==1 &
                                                                          all_diff_meth_geneIDs_with_category$`Hypermethylated in males compared to queens` ==1])
colnames(hyper_MvWQ) <- "gene_id"
write.table(hyper_MvWQ, file="diff_meth_gene_lists/GO_diff_meth_hyper_MvWQ.txt", sep="\t",quote = F, col.names = T, row.names = F)
