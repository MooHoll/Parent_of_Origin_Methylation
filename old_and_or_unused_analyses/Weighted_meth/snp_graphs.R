# Graph for SNPs 

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Methylation/Alignment_Bter")

library(ggplot2)
library(reshape2)

dat<-read.csv("snp_data.csv", header=T)
dat$hererozygous <- dat$heterozygous_0.1+dat$heterozygous_1.2
dat <- dat[,c(1,2,5,6)]
dat_subset<-reshape2::melt(dat)
head(dat_subset)

ggplot(data=dat_subset, aes(x=sample, y=value, fill=variable))+
  geom_bar(stat="identity", position="fill")+
  ylab("SNP Count")+
  ggtitle("Number of SNPs per Sample")+
  scale_fill_manual(breaks=c("filtered_snps","homozygous_alt","unique"),
                    labels = c("Total SNPs","Homozygous Alternative","Unique"),
                    values= c("dodgerblue","skyblue2","turquoise"))

  
ggplot(data=dat, aes(x=sample))+
  geom_bar(aes(y=filtered_snps), stat="identity",fill="dodgerblue",show.legend=T)+
  geom_bar(aes(y=homozygous_alt), stat="identity",fill="skyblue2",show.legend=T)+
  geom_bar(aes(y=unique), stat="identity",fill="turquoise",show.legend=T)+
  ylab("SNP Count")+
  geom_text(data = dat, aes(label = filtered_snps, y = filtered_snps -5000),
            size=5)+
  geom_text(data=dat, aes(label=unique, y=unique -5000),
            size=5)+
  #geom_text(data = dat, aes(label = homozygous_alt, y = homozygous_alt - 2500),
  #          size=5)+
  theme_bw()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20))+
  xlab("Sample")
  
## Made both graphs and put in powerpoint to put the legend of the first graph onto the second



males <- dat[1:4,]
females <- dat[5:8,]

sd(males$Total_SNPs) #144951
sd(females$Total_SNPs) #167392
