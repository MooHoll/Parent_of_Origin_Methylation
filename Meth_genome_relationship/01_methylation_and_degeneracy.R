setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/degenerate_sites")

library(readr)
library(ggplot2)
library(tidyr)
library(reshape2)

# binomial methylation calls
methylation_calls_per_sample <- read_delim("methylation_calls_per_sample.txt", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE)
head(methylation_calls_per_sample) # ~2500 CpGs methylated per sample

# genome degenreacy (one at a time as files are huge - like 3Gb each)
degenerate_sites <- read_delim("q37_degenerate_sites.txt.gz", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)
head(degenerate_sites)
degenerate_sites <- degenerate_sites[,-c(2,4)]

# Keep only C/G sites
degenerate_sites <- degenerate_sites[degenerate_sites$codon.base == "C" |
                                       degenerate_sites$codon.base == "G",]

colnames(degenerate_sites)<-c("chr","position","gene_id","codon","base","amino_acid",
                                     "codon_position","degeneracy")

# Change this depending on sample
sample <- methylation_calls_per_sample[,c(1,2,10)]
# Merge by chromosome and CpG position
all_data <- merge(sample, degenerate_sites, by=c("chr","position"))
colnames(all_data)[3] <- "sample"

all_data$sample[all_data$sample==0] <- "Unmethylated"
all_data$sample[all_data$sample==1] <- "Methylated"

# Take a look at degeneracy
count_degeneracy <- as.data.frame(table(all_data[,c(3,9)]))
count_degeneracy <- spread(count_degeneracy, key = sample, value = Freq)
count_degeneracy$proportion <- count_degeneracy$Methylated / count_degeneracy$Unmethylated

prop.test(count_degeneracy$Methylated, count_degeneracy$Unmethylated)

# Check relationship between degen and codon position
degen_codon_positon <- as.data.frame(table(all_data[,c(8,9)]))
degen_codon_positon <- spread(degen_codon_positon, key = codon_position, value = Freq)


# After running all the code above for each sample, put together for one
# more representative figure
proportion_meth_degen <- read_csv("proportion_meth_degen.csv")
head(proportion_meth_degen)
prop_meth_degen_melt <- melt(proportion_meth_degen)
prop_meth_degen_melt$Sample <- substring(prop_meth_degen_melt$Sample,1, nchar(prop_meth_degen_melt$Sample)-2)

ggplot(prop_meth_degen_melt, aes(x=variable, y =value,fill=Sample))+
  geom_boxplot()+
  xlab("Codon Degeneracy")+
  ylab("Proportion Methylated")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("Male ","Queen "),
                    values = c("midnightblue","#CC6677"))

library(multcomp)
head(prop_meth_degen_melt)
model1<-lm(value ~ Sample * variable , data=prop_meth_degen_melt)
model2<-lm(value ~ variable  , data=prop_meth_degen_melt)
anova(model1,model2) # Sig interaction
summary.lm(model1) # Everything sig
posthoc<- glht(model1, linfct = mcp(variable = 'Tukey'))
summary(posthoc) # only 0x is sig diff to the other three


# Of less interest
# ---------------------------------------------------------------

# Look at codon position
count_codon_positon <- as.data.frame(table(all_data[,c(3,8)]))
count_codon_positon <- spread(count_codon_positon, key = sample, value = Freq)
count_codon_positon$proportion <- count_codon_positon$Methylated / count_codon_positon$Unmethylated

prop.test(count_codon_positon$Methylated, count_codon_positon$Unmethylated)

# Various graphs
# ---------------------------------------------------------------

ggplot(count_codon_positon, aes(x=codon_position, y =proportion))+
  geom_bar(stat="identity")+
  xlab("Codon Position")+
  ylab("Proportion Methylated")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position ="none")

ggplot(count_codon_positon, aes(x=codon_position, y =proportion))+
  geom_bar(stat="identity")+
  xlab("Codon Position")+
  ylab("Proportion Methylated")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position ="none")

ggplot(count_degeneracy, aes(x=degeneracy, y =proportion))+
  geom_bar(stat="identity")+
  xlab("Codon Degeneracy")+
  ylab("Proportion Methylated")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position ="none")


