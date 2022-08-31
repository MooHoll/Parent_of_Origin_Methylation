## -------------------------------------------------------------------------
# Check meth levels of diff methylated genes
## ------------------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/weighted_meth")

library(readr)
library(data.table)
library(ggplot2)

# Weighted methylation and levels of all genes for each caste
weighted_meth_annotation_by_caste <- read_delim("weighted_meth_annotation_by_caste.txt", 
                                                delim = "\t", escape_double = FALSE, 
                                                trim_ws = TRUE)
head(weighted_meth_annotation_by_caste)
genes <- weighted_meth_annotation_by_caste[weighted_meth_annotation_by_caste$feature=="gene",]
genes <- genes[,c(2,6,7,8)]

genes_long <- melt(genes, id.vars = "gene_id")
head(genes_long)
colnames(genes_long) <- c("gene_id","caste","methylation")

genes_long$level <- "high"
genes_long$level[genes_long$methylation ==0] <- "none"
genes_long$level[genes_long$methylation >0 & genes_long$methylation <0.3] <- "low"
genes_long$level[genes_long$methylation >=0.3 & genes_long$methylation <0.7] <- "medium"

# Read in differentially methylated genes
MvQ_all_diff_meth_genes <- read_csv("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/coverage_files/MvQ_all_diff_meth_genes.txt")
MvW_all_diff_meth_genes <- read_csv("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/coverage_files/MvW_all_diff_meth_genes.txt")
WvQ_all_diff_meth_genes <- read_csv("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/coverage_files/WvQ_all_diff_meth_genes.txt")

MvQ_all_data <- genes_long[genes_long$gene_id %in% MvQ_all_diff_meth_genes$gene_id,]
MvW_all_data <- genes_long[genes_long$gene_id %in% MvW_all_diff_meth_genes$gene_id,]
WvQ_all_data <- genes_long[genes_long$gene_id %in% WvQ_all_diff_meth_genes$gene_id,]

ggplot(MvQ_all_data, aes(x=level, fill=caste))+
  geom_bar(stat="count", position = "dodge")+
  xlab("Methylation Level")+
  ylab("Number of Genes")+
  ggtitle("Male vs Queen Differentially Methylated Genes")+
  theme_bw()+
  scale_fill_manual(breaks = c("male","queen","worker"),labels=c("Male","Queen","Worker"),
                    values=c("#44AA99","#CC6677","#DDCC77"))+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())

ggplot(MvW_all_data, aes(x=level, fill=caste))+
  geom_bar(stat="count", position = "dodge")+
  xlab("Methylation Level")+
  ylab("Number of Genes")+
  ggtitle("Male vs Worker Differentially Methylated Genes")+
  theme_bw()+
  scale_fill_manual(breaks = c("male","queen","worker"),labels=c("Male","Queen","Worker"),
                    values=c("#44AA99","#CC6677","#DDCC77"))+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())

ggplot(WvQ_all_data, aes(x=level, fill=caste))+
  geom_bar(stat="count", position = "dodge")+
  xlab("Methylation Level")+
  ylab("Number of Genes")+
  ggtitle("Worker vs Queen Differentially Methylated Genes")+
  theme_bw()+
  scale_fill_manual(breaks = c("male","queen","worker"),labels=c("Male","Queen","Worker"),
                    values=c("#44AA99","#CC6677","#DDCC77"))+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())
