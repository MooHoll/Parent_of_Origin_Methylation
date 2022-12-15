## -------------------------------------------------------------------------
# Re-make meth over feature graph for diff levels of methylation
## -------------------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/weighted_meth")
library(readr)
library(doBy)
library(ggplot2)
library(reshape2)
library(dplyr)
library(Hmisc)
library(scales)
library(ggpubr)
library(UpSetR)
library(grid)
## -------------------------------------------------------------------------
annotation <- read_delim("weighted_meth_annotation_by_caste.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

annotation <- annotation[,-c(3,4,5)]

melted_annot <- melt(annotation, id.vars = c("feature","gene_id") )
colnames(melted_annot) <- c("Feature","ID","Caste","Weighted_Methylation")

# Remove rows where NA in one sex
melted_annot <- melted_annot[!is.na(melted_annot$Weighted_Methylation),]

# Remove mRNA and CDS not really informative
melted_annot <- melted_annot[!melted_annot$Feature == "CDS",]
melted_annot <- melted_annot[!melted_annot$Feature == "RNA",]
melted_annot <- melted_annot[!melted_annot$Feature == "mRNA",]
melted_annot <- melted_annot[!melted_annot$Feature == "transcript",]

## -------------------------------------------------------------------------
#### Define summary function (ref:http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
} 
## -------------------------------------------------------------------------
# Normal graph with all information

summary_all<-summarySE(melted_annot, measurevar = "Weighted_Methylation", 
                       groupvars = c("Feature","Caste"))

ggplot(summary_all, aes(x=Feature, y=Weighted_Methylation, fill=Caste))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=Weighted_Methylation-ci, ymax=Weighted_Methylation+ci),
                width=.2,
                position = position_dodge(.9))+
  theme_bw()+
  xlab("Genomic Feature")+
  ylab("Weighted Methylation Level")+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("male","queen","worker"),labels=c("Male","Queen","Worker"),
                    values=c("midnightblue","#CC6677","#DDCC77"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","intergenic","lnc_RNA"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","Intergenic","lnc RNA"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","lnc_RNA","intergenic"))

head(melted_annot)
median(melted_annot$Weighted_Methylation[melted_annot$Feature=="exon" &
                                           melted_annot$Caste=="worker"])


# Stats
library(multcomp)
head(melted_annot)
model1<-lm(Weighted_Methylation ~ Caste * Feature , data=melted_annot)
model2<-lm(Weighted_Methylation ~ Caste + Feature , data=melted_annot)
anova(model1,model2) # Sig interaction
summary.lm(model1) # Everything sig

# USE EAMONN'S NEW EMM THING TO SEE THE INTERACTION EFFECTS


## -------------------------------------------------------------------------
# Frequency of features being high/low methylated
head(melted_annot)
meth_low <- melted_annot[melted_annot$Weighted_Methylation <= 0.3,]
meth_medium <- melted_annot[melted_annot$Weighted_Methylation > 0.3 &
                              melted_annot$Weighted_Methylation <= 0.7,]
meth_high <- melted_annot[melted_annot$Weighted_Methylation > 0.7,]
meth_none <- melted_annot[melted_annot$Weighted_Methylation ==0,]

melted_annot$bins<-"low"
melted_annot$bins[melted_annot$Weighted_Methylation > 0.3 &
                    melted_annot$Weighted_Methylation <= 0.7] <-"medium"
melted_annot$bins[melted_annot$Weighted_Methylation > 0.7] <-"high"
melted_annot$bins[melted_annot$Weighted_Methylation ==0] <-"none"

#melted_annot$combined <- paste0(melted_annot$Feature, "_", melted_annot$Caste)
melted_meth_stuff_2 <- melted_annot

plot_data <- melted_meth_stuff_2[,c(1,3,5)]
plot_data <- aggregate(cbind(plot_data[0],counts=1), plot_data, length)


plot_data_prom_low <- subset(plot_data, bins =="low")
plot_data_prom_low$counts <- as.numeric(plot_data_prom_low$counts)

b1<- ggplot(plot_data_prom_low, aes(x=Feature, fill=Caste, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
 # xlab("Genomic Feature")+
#  ylab("Count")+
  ggtitle("Low Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("male","queen","worker"),labels=c("Male","Queen","Worker"),
                    values=c("midnightblue","#CC6677","#DDCC77"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","intergenic","lnc_RNA"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","Intergenic","lnc RNA"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","lnc_RNA","intergenic"))

plot_data_prom_med <- subset(plot_data, bins =="medium")
plot_data_prom_med$counts <- as.numeric(plot_data_prom_med$counts)

b2<- ggplot(plot_data_prom_med, aes(x=Feature, fill=Caste, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
  #xlab("Genomic Feature")+
  #ylab("Count")+
  ggtitle("Medium Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("male","queen","worker"),labels=c("Male","Queen","Worker"),
                    values=c("midnightblue","#CC6677","#DDCC77"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","intergenic","lnc_RNA"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","Intergenic","lnc RNA"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","lnc_RNA","intergenic"))

plot_data_prom_high <- subset(plot_data, bins =="high")
plot_data_prom_high$counts <- as.numeric(plot_data_prom_high$counts)

b3<- ggplot(plot_data_prom_high, aes(x=Feature, fill=Caste, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
 # xlab("")+
  ylab("Count")+
  ggtitle("High Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title.y=element_text(size=12),
        axis.title.x = element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("male","queen","worker"),labels=c("Male","Queen","Worker"),
                    values=c("midnightblue","#CC6677","#DDCC77"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","intergenic","lnc_RNA"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","Intergenic","lnc RNA"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","lnc_RNA","intergenic"))

## -------------------------------------------------------------------------
plot_data_prom <- subset(plot_data, bins =="none")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b4 <- ggplot(plot_data_prom, aes(x=Feature, fill=Caste, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
  # xlab("Genomic Feature")+
  #  ylab("Count")+
  ggtitle("No Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("male","queen","worker"),labels=c("Male","Queen","Worker"),
                    values=c("midnightblue","#CC6677","#DDCC77"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","intergenic","lnc_RNA"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","Intergenic","lnc RNA"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","lnc_RNA","intergenic"))

plot_data_prom <- subset(plot_data, bins =="high")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b3_1<- ggplot(plot_data_prom, aes(x=Feature, fill=Caste, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
  # xlab("")+
 # ylab("Count")+
  ggtitle("High Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("male","queen","worker"),labels=c("Male","Queen","Worker"),
                    values=c("midnightblue","#CC6677","#DDCC77"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","intergenic","lnc_RNA"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","Intergenic","lnc RNA"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","lnc_RNA","intergenic"))


## -------------------------------------------------------------------------
# Slightly diff figures
levels_across_feature <- ggarrange(b3_1,b2,b1,b4,
                                   ncol=2, nrow=2, common.legend = TRUE, legend="right")

annotate_figure(levels_across_feature, 
                 left = text_grob("Count", 
                                 color = "black", rot = 90, size=14),
                bottom = text_grob("Genomic Feature", 
                                   color = "black", size =12,
                                   hjust = 0.75 ))

## -------------------------------------------------------------------------
# Pull out high/low gene lists
tail(melted_annot)
genes <- melted_annot[melted_annot$Feature=="gene",]
head(genes)

# Lists of genes which are 'methylated'
general_methylated_male <- as.data.frame(unique(genes$ID[genes$Weighted_Methylation >0.005 &
                                             genes$Caste=="male"])) # 7206/11652
colnames(general_methylated_male) <- "gene_id"
write.table(general_methylated_male, file="male_all_methylated_genes.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

general_methylated_queen <- as.data.frame(unique(genes$ID[genes$Weighted_Methylation >0.005 &
                                                           genes$Caste=="queen"])) # 7754/11652
colnames(general_methylated_queen) <- "gene_id"
write.table(general_methylated_queen, file="queen_all_methylated_genes.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

general_methylated_worker <- as.data.frame(unique(genes$ID[genes$Weighted_Methylation >0.005 &
                                                            genes$Caste=="worker"])) # 7316/11652
colnames(general_methylated_worker) <- "gene_id"
write.table(general_methylated_worker, file="worker_all_methylated_genes.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

# Pull out highly methylated genes only
high_methylated_male <- as.data.frame(unique(genes$ID[genes$bins=="high" &
                                                           genes$Caste=="male"])) # 97/7206
colnames(high_methylated_male) <- "gene_id"
write.table(high_methylated_male, file="male_high_methylated_genes.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

high_methylated_queen <- as.data.frame(unique(genes$ID[genes$bins=="high" &
                                                        genes$Caste=="queen"])) # 51/7754
colnames(high_methylated_queen) <- "gene_id"
write.table(high_methylated_queen, file="queen_high_methylated_genes.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

high_methylated_worker <- as.data.frame(unique(genes$ID[genes$bins=="high" &
                                                        genes$Caste=="worker"])) # 37/7316
colnames(high_methylated_worker) <- "gene_id"
write.table(high_methylated_worker, file="worker_high_methylated_genes.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

# Quick look to see if highly methylated genes are common amongst castes/sexes
all <- rbind(high_methylated_male, high_methylated_queen, high_methylated_worker)
all <- as.data.frame(all[!duplicated(all),])
colnames(all) <- "gene_id"

all$`High methylation in males` <- 0
all$`High methylation in males`[all$gene_id %in% high_methylated_male$gene_id] <- 1

all$`High methylation in queens` <- 0
all$`High methylation in queens`[all$gene_id %in% high_methylated_queen$gene_id] <- 1

all$`High methylation in workers` <- 0
all$`High methylation in workers`[all$gene_id %in% high_methylated_worker$gene_id] <- 1

upset(all, order.by = "freq",
      text.scale = 2,
      point.size = 4,
      scale.sets = "identity",
      mainbar.y.label =NULL)
grid.text("Intersection Size",x = 0.35, y=0.60, gp=gpar(fontsize=20), rot = 90)




