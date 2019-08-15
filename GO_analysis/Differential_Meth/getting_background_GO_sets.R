# Getting base GO set for all methylated genes 

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Methylation/Diff_Meth/queens_vs_males/GO_analysis")

library(readr)

GO_annotations <- read.table("Bumble_bee_ensemble_GO_terms.txt")
colnames(GO_annotations) <- c("geneID","goID")


methylated_genes <- read_csv("gene_lists/methylated_genes_leuven_meth.txt") # 3550
colnames(methylated_genes) <- "geneID"
diff_meth_genes <- read_csv("gene_lists/diff_meth_weighted_genes.txt", col_names = F) #344
colnames(diff_meth_genes) <- "geneID"

meth_gos <- merge(GO_annotations, methylated_genes, by="geneID")
length(meth_gos[is.na(meth_gos$geneID)]) # 0
write.table(meth_gos, file = "GOs_all_methylated_genes_as_background.txt", sep="\t",
            col.names = T, row.names = F, quote = F)

diff_gos <- merge(GO_annotations, diff_meth_genes, by="geneID")
length(diff_gos[is.na(diff_gos$geneID)]) # 0
write.table(diff_gos, file = "GOs_all_diffmeth_genes_as_background.txt", sep="\t",
            col.names = T, row.names = F, quote = F)


