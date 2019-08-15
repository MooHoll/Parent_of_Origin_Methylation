# Making graphs 

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Methylation/PoO_meth/scatter_graphs")

library(readr)
library(sqldf)
library(doBy)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(FSA)
library(DescTools)
library(rcompanion)
library(multcompView)

# Data
meth_melted <- read_delim("weighted_meth/melted_weighted_meth_by_sample_type.txt", 
                                                  "\t", escape_double = FALSE, trim_ws = TRUE)
meth_melted<- meth_melted[,-7]

meth <- read_delim("weighted_meth/weighted_meth_by_sample_type.txt", 
                                           "\t", escape_double = FALSE, trim_ws = TRUE)

# Genes that are diff meth between worker alignments
PoO_meth_genes <- read_csv("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Methylation/PoO_meth/PoO_meth_genes.txt", 
                           col_names = F)

# Genes that are diff between queen and male
diff_meth_QM <- read_csv("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Methylation/Diff_Meth/queens_vs_males/list_DMGs_MSCfilter_Male_Queen.csv", 
                         col_names = FALSE)


# Remove rows where there are no values for either male/queen/workers etc.
meth_conservative <- meth_melted[complete.cases(meth_melted),] #1077 genes

# Cheeky boxplot of weighted methylation by sample group
plot_data <- subset(meth, origin == c("male","queens","workers"))

ggplot(plot_data, aes(x=origin, y=weightedMeth.mean, fill=origin))+
  geom_boxplot()+
  theme_bw()+
  xlab("Sample Type")+
  ylab("Mean Weighted Methylation Level per Gene")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.title=element_text(size = 18))+
  scale_x_discrete(labels = c('Male','Queen','Reproductive \n Worker'))+
  guides(fill=FALSE)+
  scale_fill_manual(values=c("#FF6600","#9900CC","#00CC66"))

hist(plot_data$weightedMeth.mean[plot_data$origin=="male"])
hist(plot_data$weightedMeth.mean[plot_data$origin=="queens"])
hist(plot_data$weightedMeth.mean[plot_data$origin=="workers"])

plot_data$logData <- log10(plot_data$weightedMeth.mean)

hist(plot_data$logData[plot_data$origin=="male"])
hist(plot_data$logData[plot_data$origin=="queens"])
hist(plot_data$logData[plot_data$origin=="workers"])

kruskal.test(plot_data, weightedMeth.mean ~ origin)
# Kruskal-Wallis chi-squared = 10011, df = 3, p-value < 2.2e-16

output= dunnTest( weightedMeth.mean ~ origin,
              data=plot_data,
              method="bh")
output

mean(plot_data$weightedMeth.mean[plot_data$origin=="male"])
mean(plot_data$weightedMeth.mean[plot_data$origin=="queens"])
mean(plot_data$weightedMeth.mean[plot_data$origin=="workers"])

sd(plot_data$weightedMeth.mean[plot_data$origin=="male"])
sd(plot_data$weightedMeth.mean[plot_data$origin=="queens"])
sd(plot_data$weightedMeth.mean[plot_data$origin=="workers"])

# New column to show meth difference between male/queen alignments (males-queens)
# Negative values mean higher meth in queen than male
meth_conservative$worker_allele_diff <- meth_conservative$workers_to_males -  meth_conservative$workers_to_queens

# New column to show male/queen diff
# Negative values mean higher meth in queen than male
meth_conservative$male_queen_diff <- meth_conservative$male - meth_conservative$queens
meth_conservative_subset <- subset(meth_conservative, !worker_allele_diff==0 )

qplot(x=worker_allele_diff,
      y=male_queen_diff,
      data=meth_conservative_subset,
      col=ifelse(geneID %in% PoO_meth_genes$X1, "red", "black"),
      size=I(3),
      pch = I(16),
      ylab = expression(paste("Parental Methylation Difference (Male-Queen)")),
      xlab = expression(paste("Worker Allelic Methylation Difference (Paternal Allele - Maternal Allele)")),
      geom="point")+
  scale_color_identity()+
  geom_vline(xintercept=0, colour=I("black")) +
  geom_hline(yintercept=0, colour=I("black")) +
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title = element_text(size=16))+
  geom_smooth(method=lm, data=meth_conservative_subset, colour="blue")+
  annotate("text",x=0.025,y=0.2, label="paste(italic(R) ^ 2, \" = 0.008623\")",size=10,parse=T)

# Model these data (this is non-significant!)
fit_both<-lm(male_queen_diff~worker_allele_diff, data=meth_conservative_subset, 
             na.action=na.exclude)
summary(fit_both)
anova(fit_both) 
# Adjusted R-squared:  0.008623 
# F-statistic: 1.618 on 1 and 70 DF,  p-value: 0.2076
  

# PoO meth against PoO gene exp
gene_exp <- read_csv("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Final_imprinting_model/Imprinting_graphs/stats_imprinting_model2_total counts per gene reproductive nonreproductive workers.csv")

gene_exp1 <- gene_exp[,c(2,11)]  
gene_exp_final <- gene_exp1[(!duplicated(gene_exp1$geneID)),]

workers_allele <- meth_conservative_subset[,c(1,7)]
workers_allele <- subset(workers_allele, !worker_allele_diff ==0)

all_data_meth_exp <- merge(gene_exp_final, workers_allele, by="geneID")
all_data_meth_exp <- all_data_meth_exp[(complete.cases(all_data_meth_exp)),]
all_data_meth_exp$pat_exp_porp <- 1- all_data_meth_exp$avgpropmatexpr

qplot(x=worker_allele_diff,
      y=pat_exp_porp,
      data=all_data_meth_exp,
      col=ifelse(geneID %in% PoO_meth_genes$X1, "red", "black"),
      size=I(3),
      pch = I(16),
      ylim = c(0,1),
      ylab = expression(paste("Average Proportion of Paternal Expression")),
      xlab = expression(paste("Worker Allelic Methylation Difference (Paternal Allele - Maternal Allele)")),
      geom="point")+
  scale_color_identity()+
  geom_vline(xintercept=0, colour=I("black")) +
  geom_hline(yintercept=0.6, colour=I("black")) +
  geom_hline(yintercept=0.4, colour=I("black")) +
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title = element_text(size=16))+
  geom_smooth(method=lm, data=all_data_meth_exp, colour="blue")+
  annotate("text",x=0.025,y=0.85, label="paste(italic(R) ^ 2, \" = -0.0119 \")",size=10,parse=T)


# Model these data (this is non-significant!)
fit_both1<-lm(pat_exp_porp~worker_allele_diff, data=all_data_meth_exp, 
             na.action=na.exclude)
summary(fit_both1)
anova(fit_both1) 
# Adjusted R-squared:  -0.0119 
# F-statistic: 0.2828 on 1 and 60 DF,  p-value: 0.5968

  