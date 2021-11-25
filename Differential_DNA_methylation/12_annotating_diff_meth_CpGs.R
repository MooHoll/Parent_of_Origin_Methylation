## -------------------------------------------------------------------------
# Filtering diff meth CpGs from methylkit with weighted meth of features
## -------------------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/Differential_methylation")
library(sqldf)
library(readr)
library(doBy)
library(ggplot2)
library(dplyr)
library(data.table)
library(sqldf)
library(grid)

annotation <- read_delim("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/ref_Bter_1.0_top_level_numbered_exons.txt", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)

diff_meth_sites_MvW <- read_csv("male_vs_worker_diff_meth_CpGs.csv")
diff_meth_sites_MvW <- diff_meth_sites_MvW[,c(2,3,8)]
colnames(diff_meth_sites_MvW) <- c("chr","cpg_position","methylation_diff")

diff_meth_sites_MvQ <- read_csv("male_vs_queen_diff_meth_CpGs.csv")
diff_meth_sites_MvQ <- diff_meth_sites_MvQ[,c(2,3,8)]
colnames(diff_meth_sites_MvQ) <- c("chr","cpg_position","methylation_diff")

diff_meth_sites_WvQ <- read_csv("worker_vs_queen_diff_meth_CpGs.csv")
diff_meth_sites_WvQ <- diff_meth_sites_WvQ[,c(2,3,8)]
colnames(diff_meth_sites_WvQ) <- c("chr","cpg_position","methylation_diff")

## -------------------------------------------------------------------------

# Add column for hypermeth sex
diff_meth_sites_MvW$hypermethylated <- "male"
diff_meth_sites_MvW$hypermethylated[diff_meth_sites_MvW$methylation_diff >0] <- "worker"

diff_meth_sites_MvQ$hypermethylated <- "male"
diff_meth_sites_MvQ$hypermethylated[diff_meth_sites_MvQ$methylation_diff >0] <- "queen"

diff_meth_sites_WvQ$hypermethylated <- "worker"
diff_meth_sites_WvQ$hypermethylated[diff_meth_sites_WvQ$methylation_diff >0] <- "queen"

# How many hypermeth in each sex?
nrow(diff_meth_sites_MvW[diff_meth_sites_MvW$hypermethylated=="male",]) #501
nrow(diff_meth_sites_MvW[diff_meth_sites_MvW$hypermethylated=="worker",]) #510

nrow(diff_meth_sites_MvQ[diff_meth_sites_MvQ$hypermethylated=="male",]) #369
nrow(diff_meth_sites_MvQ[diff_meth_sites_MvQ$hypermethylated=="queen",]) #455

nrow(diff_meth_sites_WvQ[diff_meth_sites_WvQ$hypermethylated=="worker",]) #84
nrow(diff_meth_sites_WvQ[diff_meth_sites_WvQ$hypermethylated=="queen",]) #72

# Goodness of fit
observed = c(84, 72)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions

chisq.test(x = observed,
           p = expected)
# male vs worker: X-squared = 0.080119, df = 1, p-value = 0.7771
# male vs queen: X-squared = 8.9757, df = 1, p-value = 0.002736
# worker vs queen: X-squared = 0.92308, df = 1, p-value = 0.3367

diff_meth_sites_MvW$comparison <- "MvW"
diff_meth_sites_MvQ$comparison <- "MvQ"
diff_meth_sites_WvQ$comparison <- "WvQ" 

one_df <- rbind(diff_meth_sites_MvW,diff_meth_sites_MvQ,diff_meth_sites_WvQ)

## -------------------------------------------------------------------------
# Annotate the differentially methylated CpGs with genomic features

output <- sqldf("SELECT sample.chr,
                    sample.cpg_position,
                    sample.methylation_diff,
                    sample.hypermethylated,
                    sample.comparison,
                    annot.chr,
                    annot.feature,
                    annot.start,
                    annot.end,
                    annot.gene_id,
                    annot.number
                    FROM one_df AS sample
                    LEFT JOIN annotation AS annot
                    ON sample.chr = annot.chr
                    AND (sample.cpg_position >= annot.start AND sample.cpg_position <= annot.end)")

output <- output[,-1]

## -------------------------------------------------------------------------
# Where are these CpGs
## -------------------------------------------------------------------------

not_in_feature <- output[is.na(output$gene_id),] #53/8088 not in feature
output$feature <- as.factor(output$feature)

# Remove non-informative locations
output <- output[!output$feature == "CDS",]
output <- output[!output$feature == "gene",]
output <- output[!output$feature == "mRNA",]
output <- output[!output$feature == "RNA",]
output <- output[!output$feature == "transcript",]
output <- output[!is.na(output$feature),]
# Note: 2579 annotations from 1991 positions as some positions fall over multiple annotations

ggplot(output, aes(x=feature, fill=hypermethylated))+
  geom_bar()+
  guides()+
  xlab("Feature")+
  ylab("Number of Significant CpGs")+
  theme_bw()+
  facet_grid(~comparison)+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("male","queen","worker"),labels=c("Male","Queen","Worker"),
                    values=c("#44AA99","#CC6677","#DDCC77"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","intergenic","lnc_RNA"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","Intergenic","lnc RNA"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","lnc_RNA","intergenic"))


## -------------------------------------------------------------------------
# average number diff CpGs by feature
## -------------------------------------------------------------------------

head(output)
output_count <- summaryBy(cpg_position ~ feature + gene_id + comparison , data=output, FUN=length)
output_count <- output_count[!output_count$feature=="intergenic",]

ggplot(output_count, aes(x=feature, y=cpg_position.length))+
  geom_boxplot()+
  xlab("Feature")+
  ylab("Number of Significant CpGs")+
  theme_bw()+
  facet_grid(~comparison)+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 14))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","intergenic","lnc_RNA"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","Intergenic","lnc RNA"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","lnc_RNA","intergenic"))

# Here we are saying a min of 2 diff meth cpg per feature to count
exon_data <- output[output$feature=="exon",]
number_diff_cpgs_per_exon <- dplyr::count(exon_data, gene_id)
range(number_diff_cpgs_per_exon$n) #1-23, mean = 2.2, median = 1
hist(number_diff_cpgs_per_exon$n)
nrow(number_diff_cpgs_per_exon[number_diff_cpgs_per_exon$n >=2,]) #354
exons_with_2_cpgs <- subset(number_diff_cpgs_per_exon, n >=2)

# Write out the exon info for later scatter
for_scatter <- exon_data[exon_data$gene_id %in% exons_with_2_cpgs$gene_id,]
for_scatter <- for_scatter[,c(3,4,7,9,10)]
write.table(for_scatter, file="diff_meth_exons_list.txt", quote = F,
            col.names = T, row.names = F)

# Which number exons are these CpGs in?
head(exon_data)

# Are these CpGs localised to the first few exons?
ggplot(exon_data, aes(x=number))+
  geom_bar(stat = "count")+
  xlab("Exon Number")+
  ylab("Number of Significant CpGs")+
  theme_bw()+
  facet_grid(~comparison)+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())

## -------------------------------------------------------------------------
# After have lists, filter on weighted meth difference of exon/gene
## -------------------------------------------------------------------------

head(exons_with_2_cpgs)
exon_data <- output[(output$gene_id %in% exons_with_2_cpgs$gene_id &
                         output$feature == "exon"),]

# Read in the weighted methylation levels
weighted_meth <- read_delim("weighted_meth/weighted_meth_annotation_by_caste.txt", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)
head(weighted_meth)
weighted_meth$MvW_diff <- weighted_meth$worker - weighted_meth$male
weighted_meth$MvQ_diff <- weighted_meth$queen - weighted_meth$male
weighted_meth$WvQ_diff <- weighted_meth$queen - weighted_meth$worker

# remove rows where no data for male or female
weighted_meth <- weighted_meth[!is.na(weighted_meth$male),]
weighted_meth <- weighted_meth[!is.na(weighted_meth$queen),]
weighted_meth <- weighted_meth[!is.na(weighted_meth$worker),]

# Keep only features where the meth level is greatre than background 0.05 in at least one sex
weighted_meth <- weighted_meth[(weighted_meth$male > 0.005 | weighted_meth$queen > 0.005 |
                                  weighted_meth$worker > 0.005),]

# Keep genes that have at least two diff CpGs in exons and 15% feature level difference overall
MvW_exons <- exon_data[exon_data$comparison=="MvW",]
weighted_meth_exons_MvW <- weighted_meth[(weighted_meth$gene_id %in% MvW_exons$gene_id &
                                       weighted_meth$feature == "exon"),] 
weighted_meth_exons_MvW <- weighted_meth_exons_MvW[(weighted_meth_exons_MvW$MvW_diff > 0.15 |
                                                      weighted_meth_exons_MvW$MvW_diff < -0.15),]
write.table(weighted_meth_exons_MvW,file="weighted_meth_exons_MvW.txt",sep="\t",quote = F,
            col.names = T, row.names = F)
length(unique(weighted_meth_exons_MvW$gene_id)) # 155 genes
genes <- as.data.frame(unique(weighted_meth_exons_MvW$gene_id))
colnames(genes)<-"gene_id"
write.table(genes, file="MvW_all_diff_meth_genes.txt",sep="\t",quote = F,
            col.names = T, row.names = F)


MvQ_exons <- exon_data[exon_data$comparison=="MvQ",]
weighted_meth_exons_MvQ <- weighted_meth[(weighted_meth$gene_id %in% MvQ_exons$gene_id &
                                            weighted_meth$feature == "exon"),] 
weighted_meth_exons_MvQ <- weighted_meth_exons_MvQ[(weighted_meth_exons_MvQ$MvQ_diff > 0.15 |
                                                      weighted_meth_exons_MvQ$MvQ_diff < -0.15),]
write.table(weighted_meth_exons_MvQ,file="weighted_meth_exons_MvQ.txt",sep="\t",quote = F,
            col.names = T, row.names = F)
length(unique(weighted_meth_exons_MvQ$gene_id)) # 165 genes
genes <- as.data.frame(unique(weighted_meth_exons_MvQ$gene_id))
colnames(genes)<-"gene_id"
write.table(genes, file="MvQ_all_diff_meth_genes.txt",sep="\t",quote = F,
            col.names = T, row.names = F)

WvQ_exons <- exon_data[exon_data$comparison=="WvQ",]
weighted_meth_exons_WvQ <- weighted_meth[(weighted_meth$gene_id %in% WvQ_exons$gene_id &
                                            weighted_meth$feature == "exon"),] 
weighted_meth_exons_WvQ <- weighted_meth_exons_WvQ[(weighted_meth_exons_WvQ$WvQ_diff > 0.15 |
                                                      weighted_meth_exons_WvQ$WvQ_diff < -0.15),]
write.table(weighted_meth_exons_WvQ,file="weighted_meth_exons_WvQ.txt",sep="\t",quote = F,
            col.names = T, row.names = F)
length(unique(weighted_meth_exons_WvQ$gene_id)) # 37 genes
genes <- as.data.frame(unique(weighted_meth_exons_WvQ$gene_id))
colnames(genes)<-"gene_id"
write.table(genes, file="WvQ_all_diff_meth_genes.txt",sep="\t",quote = F,
            col.names = T, row.names = F)


# Which sex are these genes hypermethylated in
weighted_meth_exons_MvW$hypermethylated <- "male"
weighted_meth_exons_MvW$hypermethylated[weighted_meth_exons_MvW$worker > weighted_meth_exons_MvW$male] <- "worker"

male_hyper_exon_genes <- unique(weighted_meth_exons_MvW$gene_id[weighted_meth_exons_MvW$hypermethylated=="male"])
length(male_hyper_exon_genes) # 95
genes <- as.data.frame(male_hyper_exon_genes)
colnames(genes)<-"gene_id"
write.table(genes, file="MvW_male_hypermeth_genes.txt",sep="\t",quote = F,
            col.names = T, row.names = F)

worker_hyper_exon_genes <- unique(weighted_meth_exons_MvW$gene_id[weighted_meth_exons_MvW$hypermethylated=="worker"])
length(worker_hyper_exon_genes) # 107
genes <- as.data.frame(worker_hyper_exon_genes)
colnames(genes)<-"gene_id"
write.table(genes, file="MvW_worker_hypermeth_genes.txt",sep="\t",quote = F,
            col.names = T, row.names = F)

both <- Reduce(intersect, list(male_hyper_exon_genes,worker_hyper_exon_genes)) # 47
common_exon_MvW <- weighted_meth_exons_MvW[weighted_meth_exons_MvW$gene_id %in% both,]


weighted_meth_exons_MvQ$hypermethylated <- "male"
weighted_meth_exons_MvQ$hypermethylated[weighted_meth_exons_MvQ$queen > weighted_meth_exons_MvQ$male] <- "queen"

male_hyper_exon_genes_2 <- unique(weighted_meth_exons_MvQ$gene_id[weighted_meth_exons_MvQ$hypermethylated=="male"])
length(male_hyper_exon_genes_2) # 94
genes <- as.data.frame(male_hyper_exon_genes_2)
colnames(genes)<-"gene_id"
write.table(genes, file="MvQ_male_hypermeth_genes.txt",sep="\t",quote = F,
            col.names = T, row.names = F)

queen_hyper_exon_genes <- unique(weighted_meth_exons_MvQ$gene_id[weighted_meth_exons_MvQ$hypermethylated=="queen"])
length(queen_hyper_exon_genes) # 120
genes <- as.data.frame(queen_hyper_exon_genes)
colnames(genes)<-"gene_id"
write.table(genes, file="MvQ_queen_hypermeth_genes.txt",sep="\t",quote = F,
            col.names = T, row.names = F)

both <- Reduce(intersect, list(male_hyper_exon_genes_2,queen_hyper_exon_genes)) # 49
common_exon_MvQ <- weighted_meth_exons_MvQ[weighted_meth_exons_MvQ$gene_id %in% both,]


weighted_meth_exons_WvQ$hypermethylated <- "worker"
weighted_meth_exons_WvQ$hypermethylated[weighted_meth_exons_WvQ$queen > weighted_meth_exons_WvQ$worker] <- "queen"

worker_hyper_exon_genes_2 <- unique(weighted_meth_exons_WvQ$gene_id[weighted_meth_exons_WvQ$hypermethylated=="worker"])
length(worker_hyper_exon_genes_2) # 27
genes <- as.data.frame(worker_hyper_exon_genes_2)
colnames(genes)<-"gene_id"
write.table(genes, file="WvQ_worker_hypermeth_genes.txt",sep="\t",quote = F,
            col.names = T, row.names = F)

queen_hyper_exon_genes_2 <- unique(weighted_meth_exons_WvQ$gene_id[weighted_meth_exons_WvQ$hypermethylated=="queen"])
length(queen_hyper_exon_genes_2) # 19
genes <- as.data.frame(queen_hyper_exon_genes_2)
colnames(genes)<-"gene_id"
write.table(genes, file="WvQ_queen_hypermeth_genes.txt",sep="\t",quote = F,
            col.names = T, row.names = F)
            
both <- Reduce(intersect, list(worker_hyper_exon_genes_2,queen_hyper_exon_genes_2)) # 9
common_exon_WvQ <- weighted_meth_exons_WvQ[weighted_meth_exons_WvQ$gene_id %in% both,]

# Make a plot to make this easier to interpret
weighted_meth_exons_MvW$comparison <- "MvW"
weighted_meth_exons_MvW$hypermethylated[weighted_meth_exons_MvW$gene_id %in% common_exon_MvW$gene_id] <- "both"
weighted_meth_exons_MvQ$comparison <- "MvQ"
weighted_meth_exons_MvQ$hypermethylated[weighted_meth_exons_MvQ$gene_id %in% common_exon_MvQ$gene_id] <- "both"
weighted_meth_exons_WvQ$comparison <- "WvQ"
weighted_meth_exons_WvQ$hypermethylated[weighted_meth_exons_WvQ$gene_id %in% common_exon_WvQ$gene_id] <- "both"

all_for_plot <- rbind(weighted_meth_exons_MvW,weighted_meth_exons_MvQ,weighted_meth_exons_WvQ)
all_for_plot <- all_for_plot[,c(2,12,13)]
all_for_plot <- all_for_plot[!duplicated(all_for_plot),]

ggplot(all_for_plot, aes(x=comparison, fill=hypermethylated))+
  geom_bar()+
  guides()+
  xlab("Comparison")+
  ylab("Number of Genes")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("male","queen","worker","both"),labels=c("Male","Queen","Worker","Both"),
                    values=c("#44AA99","#CC6677","#DDCC77","grey64"))+
  scale_x_discrete(breaks=c("MvW","MvQ","WvQ"),
                   labels=c("Male vs Worker", "Male vs Queen","Worker vs Queen"))

# Oooooh a better way is to do an upset plot!!! Fuck yes! 
library(UpSetR)

hypermethylated_MvW_males <- as.data.frame(weighted_meth_exons_MvW$gene_id[weighted_meth_exons_MvW$hypermethylated=="male"])
colnames(hypermethylated_MvW_males)<-"gene_id"
hypermethylated_MvW_workers <- as.data.frame(weighted_meth_exons_MvW$gene_id[weighted_meth_exons_MvW$hypermethylated=="worker"])
colnames(hypermethylated_MvW_workers)<-"gene_id"

hypermethylated_MvQ_males <- as.data.frame(weighted_meth_exons_MvQ$gene_id[weighted_meth_exons_MvQ$hypermethylated=="male"])
colnames(hypermethylated_MvQ_males)<-"gene_id"
hypermethylated_MvQ_queens <- as.data.frame(weighted_meth_exons_MvQ$gene_id[weighted_meth_exons_MvQ$hypermethylated=="queen"])
colnames(hypermethylated_MvQ_queens)<-"gene_id"

hypermethylated_WvQ_workers <- as.data.frame(weighted_meth_exons_WvQ$gene_id[weighted_meth_exons_WvQ$hypermethylated=="worker"])
colnames(hypermethylated_WvQ_workers)<-"gene_id"
hypermethylated_WvQ_queens <- as.data.frame(weighted_meth_exons_WvQ$gene_id[weighted_meth_exons_WvQ$hypermethylated=="queen"])
colnames(hypermethylated_WvQ_queens)<-"gene_id"

all <- rbind(hypermethylated_MvW_males,hypermethylated_MvW_workers,hypermethylated_MvQ_males,
             hypermethylated_MvQ_queens,hypermethylated_WvQ_workers,hypermethylated_WvQ_queens)
all <- as.data.frame(all[!duplicated(all),]) #178 genes which change across castes and sexes
colnames(all) <- "gene_id"

all$`Hypermethylated in males compared to workers` <- 0
all$`Hypermethylated in males compared to workers`[all$gene_id %in% hypermethylated_MvW_males$gene_id] <- 1
all$`Hypermethylated in workers compared to males` <- 0
all$`Hypermethylated in workers compared to males`[all$gene_id %in% hypermethylated_MvW_workers$gene_id] <- 1

all$`Hypermethylated in males compared to queens` <- 0
all$`Hypermethylated in males compared to queens`[all$gene_id %in% hypermethylated_MvQ_males$gene_id] <- 1
all$`Hypermethylated in queens compared to males` <- 0
all$`Hypermethylated in queens compared to males`[all$gene_id %in% hypermethylated_MvQ_queens$gene_id] <- 1

all$`Hypermethylated in workers compared to queens` <- 0
all$`Hypermethylated in workers compared to queens`[all$gene_id %in% hypermethylated_WvQ_workers$gene_id] <- 1
all$`Hypermethylated in queens compared to workers` <- 0
all$`Hypermethylated in queens compared to workers`[all$gene_id %in% hypermethylated_WvQ_queens$gene_id] <- 1

upset(all, nsets =6, order.by = "freq",
      text.scale = 2,
      point.size = 4,
      scale.sets = "identity",
      mainbar.y.label =NULL)
grid.text("Intersection Size",x = 0.545, y=0.60, gp=gpar(fontsize=20), rot = 90)

## -------------------------------------------------------------------------
# Write out all the gene lists for later use
## -------------------------------------------------------------------------

### UP TO HERE!
head(all)
write.table(all, file="all_diff_meth_geneIDs_with_category.txt", sep="\t", quote = F,
            col.names = T, row.names = F)











