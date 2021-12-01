## -------------------------------------------------------------------------
# Overlap ASM in mined data with PoO meth from Wu et al.
## -------------------------------------------------------------------------

library(readr)
library(tidyr)
library(data.table)
library(dplyr)
library(ggplot2)
library(doBy)

setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/ASM_honeybee/asm_in_genes")

# Read in asm calls (91 files)
file.list = list.files(("./"),pattern="*genes.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)

# Add the sample name as the dataframe name
base <- tools::file_path_sans_ext(basename(file.list))
base <- gsub("_asm_in_genes","",base)
names(samples) <- base

# Read in the parent-of-origin genes
parent_of_origin_genes_HB <- read_delim("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/ASM_honeybee/parent_of_origin_genes_HB.txt", 
                                        delim = "\t", escape_double = FALSE, 
                                        trim_ws = TRUE)

for(i in seq_along(samples)){
  # Fix the gene name columns
  samples[[i]] <- separate(data = samples[[i]], col = gene_id, into = c("beebase", "LOC"), sep = "\\,")
  samples[[i]]$LOC <- ifelse(is.na(samples[[i]]$LOC), samples[[i]]$beebase, samples[[i]]$LOC)
  samples[[i]]$beebase <- ifelse(samples[[i]]$beebase %like% "GeneID", NA, samples[[i]]$beebase)
  samples[[i]]$beebase <- sub("BEEBASE:","", samples[[i]]$beebase)
  samples[[i]]$LOC <- sub("GeneID=","", samples[[i]]$LOC)
  colnames(samples[[i]])[8] <- "gene_id"
  # Add a column for sample name
  samples[[i]]$sample_name <- paste0(names(samples[i]))
}
## -------------------------------------------------------------------------

head(samples[[1]])
for_plotting <- as.data.frame(bind_rows(samples))
for_plotting$sample_name <- gsub("trim_","",for_plotting$sample_name)

# Add tissue type/study metadata and remove samples which didn't pass quality checks
sample_metadata <- read_delim("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/ASM_honeybee/sample_metadata.txt", 
                              delim = "\t", escape_double = FALSE, 
                              col_names = FALSE, trim_ws = TRUE)
colnames(sample_metadata) <- c("sample_name","tissue")
for_plotting <- merge(sample_metadata, for_plotting, by="sample_name") # 58 datasets including males
for_plotting <- for_plotting[,c(1,2,10)]
for_plotting <- for_plotting[!duplicated(for_plotting),]

# Have a look at the number of ASM regions in genes
for_boxplot <- summaryBy(gene_id ~ tissue + sample_name, data=for_plotting, FUN=length) 
for_boxplot_female <- for_boxplot[!for_boxplot$tissue %like% "male",]

ggplot(for_boxplot_female, aes(x=tissue, y=gene_id.length, fill=tissue))+
  geom_boxplot(outlier.color ="red")+
  geom_jitter(size=3)+
  xlab("Sample Type")+
  ylab("Number of Genes with ASM")+
  theme_bw()+
  theme(axis.text.y=element_text(size=20),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title=element_text(size=26),
        legend.position = "none")

for_boxplot_male <- for_boxplot[for_boxplot$tissue %like% "male",]
ggplot(for_boxplot_male, aes(x=tissue, y=gene_id.length, fill=tissue))+
  geom_boxplot(outlier.color ="red")+
  geom_jitter(size=3)+
  xlab("Sample Type")+
  ylab("Number of Genes with ASM")+
  theme_bw()+
  theme(axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=20),
        axis.title=element_text(size=26),
        legend.position = "none")

## -------------------------------------------------------------------------

# Merge with honeybee samples to keep only ASM regions which also show parent-of-origin meth
samples_PoP <- list()
for(i in seq_along(samples)){
  samples_PoP[[i]] <- merge(samples[[i]], parent_of_origin_genes_HB, by="gene_id")
  samples_PoP[[i]] <- samples_PoP[[i]][,c(1,9,10,11)]
  samples_PoP[[i]] <- samples_PoP[[i]][!duplicated(samples_PoP[[i]]$gene_id),]
}

for_plotting <- as.data.frame(bind_rows(samples_PoP))
for_plotting$sample_name <- gsub("trim_","",for_plotting$sample_name)
for_plotting <- merge(sample_metadata, for_plotting, by="sample_name")
for_plotting <- for_plotting[,-4]
length(unique(for_plotting$gene_id)) # 134 /166 genes (132 for just female samples)

# Look at the female data only for now
for_plotting_female <- for_plotting[!for_plotting$tissue %like% "male",]

asm_number_in_PoO_genes <- as.data.frame(table(for_plotting_female$gene_id))
hist(asm_number_in_PoO_genes$Freq, breaks=33, xlab = "Number of data sets",
     ylab = "Number of genes", main="", cex.lab=1.8, cex.axis=1.8)

for_boxplot <- summaryBy(gene_id ~ tissue + sample_name, data=for_plotting_female, FUN=length) 
for_boxplot$total_PoO <- 166
for_boxplot$proportion <- for_boxplot$gene_id.length / for_boxplot$total_PoO
plot(for_boxplot$proportion) # Two distributions, one male one female

ggplot(for_boxplot, aes(x=tissue, y=proportion, fill=tissue))+
  geom_boxplot(outlier.color ="red")+
  geom_jitter(size=3)+
  xlab("Sample Type")+
  ylab("Proportion of PoO-methylation genes with ASM")+
  theme_bw()+
  theme(axis.text.y=element_text(size=20),
        axis.text.x=element_text(angle=45,hjust=1,size=17),
        axis.title=element_text(size=18),
        legend.position = "none")
range(for_boxplot$proportion)

# Have a look by the parent of origin
head(for_plotting_female)
head(asm_number_in_PoO_genes)

core_set <- subset(asm_number_in_PoO_genes, Freq==33)
colnames(core_set) <- c("gene_id","number")
core_set <- merge(core_set, parent_of_origin_genes_HB, by="gene_id")
table(core_set$parent) # Equal numbers maternal and paternal, 6 each

# Make a table of the 12 core genes for the paper
all_data <- as.data.frame(bind_rows(samples))
head(all_data)
all_data_core <- all_data[all_data$gene_id %in% core_set$gene_id,]
all_data_core <- all_data_core[all_data_core$sample_name %in% for_plotting_female$sample_name,]
mean_ASM_core_12 <- summaryBy(asm_score ~ gene_id, data=all_data_core, FUN=mean) 

second_set <- subset(asm_number_in_PoO_genes, Freq >= 16) #87
colnames(second_set) <- c("gene_id","number")
second_set <- merge(second_set, parent_of_origin_genes_HB, by="gene_id")
table(second_set$parent) # Equal numbers 46 mat and 45 pat
head(all_data)
all_data_female <- all_data[all_data$sample_name %in% for_plotting_female$sample_name,]
length(unique(all_data_female$LOC)) #3448
# Get list for GO background
GO_all_ASM <- as.data.frame(unique(all_data_female$LOC))
colnames(GO_all_ASM) <- "gene_id"
write.table(GO_all_ASM, file="all_ASM_background_female.txt", sep='\t', quote = F, col.names = T, row.names = F)

all_data_female <- all_data_female[,c(4,8,9)]
all_data_female <- merge(second_set, all_data_female, by ="gene_id")
mean_ASM_50percent <- summaryBy(asm_score ~ gene_id + LOC + number + parent, data=all_data_female, FUN=mean) 
write.table(mean_ASM_50percent, file="PoO_in_50percent_datasets.txt", sep='\t', quote = F, col.names = T, row.names = F)

## -------------------------------------------------------------------------
# Start fresh figure out what's going on with the male genes
for_males <- as.data.frame(bind_rows(samples))
for_males$sample_name <- gsub("trim_","",for_males$sample_name)

# Add tissue type/study metadata and remove samples which didn't pass quality checks
sample_metadata <- read_delim("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/ASM_honeybee/sample_metadata.txt", 
                              delim = "\t", escape_double = FALSE, 
                              col_names = FALSE, trim_ws = TRUE)
colnames(sample_metadata) <- c("sample_name","tissue")
for_males <- merge(sample_metadata, for_males, by="sample_name") # 58 datasets including males
for_males <- for_males[,c(1,2,10)]
for_males <- for_males[!duplicated(for_males),]
for_males <- for_males[for_males$tissue %like% "Adult male: thorax",]
length(unique(for_males$gene_id)) # 1838

samples_PoP <- list()
for(i in seq_along(samples)){
  samples_PoP[[i]] <- merge(samples[[i]], parent_of_origin_genes_HB, by="gene_id")
  samples_PoP[[i]] <- samples_PoP[[i]][,c(1,9,10,11)]
  samples_PoP[[i]] <- samples_PoP[[i]][!duplicated(samples_PoP[[i]]$gene_id),]
}

for_plotting <- as.data.frame(bind_rows(samples_PoP))
for_plotting$sample_name <- gsub("trim_","",for_plotting$sample_name)
for_plotting <- merge(sample_metadata, for_plotting, by="sample_name")
for_plotting <- for_plotting[,-4]
for_plotting <- for_plotting[for_plotting$tissue %like% "Adult male: thorax",]
length(unique(for_plotting$gene_id)) # 79 /166 genes
for_plotting <- for_plotting[!duplicated(for_plotting),]

for_boxplot <- summaryBy(gene_id ~ tissue + sample_name, data=for_plotting, FUN=length) 
for_boxplot$total_PoO <- 166
for_boxplot$proportion <- for_boxplot$gene_id.length / for_boxplot$total_PoO
plot(for_boxplot$proportion) 

ggplot(for_boxplot, aes(x=tissue, y=proportion, fill=tissue))+
  geom_boxplot(outlier.color ="red")+
  geom_jitter(size=3)+
  xlab("Sample Type")+
  ylab("Proportion of PoO-methylation genes with ASM")+
  theme_bw()+
  theme(axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=18),
        axis.title=element_text(size=18),
        legend.position = "none")
range(for_boxplot$proportion)

# Have a look by the parent of origin
head(for_plotting)
asm_number_in_PoO_genes <- as.data.frame(table(for_plotting$gene_id)) # 1 gene in all 10 datasest

core_set <- subset(asm_number_in_PoO_genes, Freq==33)
colnames(core_set) <- c("gene_id","number")
core_set <- merge(core_set, parent_of_origin_genes_HB, by="gene_id")
table(core_set$parent) # Equal numbers maternal and paternal, 6 each

# Make a set of genes present in 50% of the datasets
second_set <- subset(asm_number_in_PoO_genes, Freq >= 5) #12
colnames(second_set) <- c("gene_id","number")
second_set <- merge(second_set, parent_of_origin_genes_HB, by="gene_id")
table(second_set$parent) # Equal numbers 7 mat and 5 pat

all_data <- as.data.frame(bind_rows(samples))
head(all_data)

all_data_male <- all_data[all_data$sample_name %in% for_plotting$sample_name,]
length(unique(all_data_male$LOC)) #2010
# Get list for GO background
GO_all_ASM <- as.data.frame(unique(all_data_male$LOC))
colnames(GO_all_ASM) <- "gene_id"
write.table(GO_all_ASM, file="all_ASM_background_male.txt", sep='\t', quote = F, col.names = T, row.names = F)

all_data_male <- all_data_male[,c(4,8,9)]
all_data_male <- merge(second_set, all_data_male, by ="gene_id")
mean_ASM_50percent <- summaryBy(asm_score ~ gene_id + LOC + number + parent, data=all_data_male, FUN=mean) 
write.table(mean_ASM_50percent, file="PoO_in_50percent_datasets_male.txt", sep='\t', quote = F, col.names = T, row.names = F)



