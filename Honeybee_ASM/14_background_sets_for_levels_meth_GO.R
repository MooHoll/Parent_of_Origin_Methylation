#---------------------------------------------------
# Background gene sets for levels of meth GO analysis
#---------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/ASM_honeybee/GO_analysis")
library(readr)

Apis_mellifera_HGD_go_annotation <- read_delim("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/ASM_honeybee/GO_analysis/Apis_mellifera_HGD_go_annotation.gaf", 
                                               delim = "\t", escape_double = FALSE, 
                                               col_names = FALSE, trim_ws = TRUE)

Apis_mellifera_HGD_go_annotation <- Apis_mellifera_HGD_go_annotation[,c(3,5)]
head(Apis_mellifera_HGD_go_annotation)

colnames(Apis_mellifera_HGD_go_annotation) <- c("gene_id","GO_term")

#---------------------------------------------------
# ASM background list
#---------------------------------------------------

# ASM genes background sets
all_ASM_background_female <- read_csv("all_ASM_background_female.txt")
# 3448 genes

all_ASM_background_female_GO <- merge(all_ASM_background_female, Apis_mellifera_HGD_go_annotation, by = "gene_id")
length(unique(all_ASM_background_female_GO$gene_id)) #3181/3448

write.table(all_ASM_background_female_GO, file="GO_all_female_ASM_genes_from_HGD.txt", sep="\t",quote = F, col.names = T, row.names = F)


# For male background
all_ASM_background_male <- read_csv("all_ASM_background_male.txt")
# 2010

all_ASM_background_male_GO <- merge(all_ASM_background_male, Apis_mellifera_HGD_go_annotation, by = "gene_id")
length(unique(all_ASM_background_male_GO$gene_id)) #1876/2010

write.table(all_ASM_background_male_GO, file="GO_all_male_ASM_genes_from_HGD.txt", sep="\t",quote = F, col.names = T, row.names = F)

#---------------------------------------------------
# Have a quick look which terms are actually associated with these 12 genes lists
core_12_list <- read_csv("core_12_list.txt",  col_names = FALSE)
colnames(core_12_list) <- "gene_id"
core_12_GO <- merge(core_12_list, Apis_mellifera_HGD_go_annotation, by = "gene_id")
length(unique(core_12_GO$gene_id)) # 12/12
write.table(core_12_GO, file="core_12_with_GO.txt", sep="\t",quote = F, col.names = T, row.names = F)

male_12_list <- read_csv("male_12.txt",  col_names = FALSE)
colnames(male_12_list) <- "gene_id"
male_12_GO <- merge(male_12_list, Apis_mellifera_HGD_go_annotation, by = "gene_id")
length(unique(male_12_GO$gene_id)) # 12/12
write.table(male_12_GO, file="male_12_with_GO.txt", sep="\t",quote = F, col.names = T, row.names = F)

# Wack these lists into revigo
core_12_all_GOs <- read_csv("core_12_all_GOs.csv")
core_12_all_GOs <- core_12_all_GOs[,c(1,2)]
colnames(core_12_all_GOs)[1]<-"GO_term"
together_core_12 <- merge(core_12_GO, core_12_all_GOs, by="GO_term")
write.table(together_core_12, file="core_12_with_GO.txt", sep="\t",quote = F, col.names = T, row.names = F)

male_12_all_GOs <- read_csv("male_12_all_GOs.csv")
male_12_all_GOs <- male_12_all_GOs[,c(1,2)]
colnames(male_12_all_GOs)[1]<-"GO_term"
together_male_12 <- merge(male_12_GO, male_12_all_GOs, by="GO_term")
write.table(together_male_12, file="male_12_with_GO.txt", sep="\t",quote = F, col.names = T, row.names = F)
