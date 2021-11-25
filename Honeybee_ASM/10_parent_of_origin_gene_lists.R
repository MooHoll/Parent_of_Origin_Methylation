# Wu et al. gene selection

# Take all unique genes from both repro and sterile from both genetic blocks

setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/ASM_honeybee")

library(readxl)

Wu_PoO_genes <- read_excel("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/ASM_honeybee/Wu_PoO_genes.xlsx")

maternal_only <- Wu_PoO_genes[Wu_PoO_genes$Bias_Class =="Maternal",]
paternal_only <- Wu_PoO_genes[Wu_PoO_genes$Bias_Class =="Paternal",]

length(unique(maternal_only$Gene_ID)) #75
length(unique(paternal_only$Gene_ID)) #91

mat_genes <- as.data.frame(unique(maternal_only$Gene_ID))
colnames(mat_genes) <- "gene_id"

pat_genes <- as.data.frame(unique(paternal_only$Gene_ID))
colnames(pat_genes) <- "gene_id"

look <- merge(mat_genes, pat_genes, by="gene_id") # 0 in common, phew.

mat_genes$parent <- "maternal"
pat_genes$parent <- "paternal"

all <- rbind(mat_genes, pat_genes)
write.table(all, file="parent_of_origin_genes_HB.txt", sep="\t", quote = F,
            col.names = T, row.names = F)
