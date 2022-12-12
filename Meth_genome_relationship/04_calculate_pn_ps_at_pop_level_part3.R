## -------------------------------------------------------------------------
# Calculate pn/ps
## ------------------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/labelled_snps")

library(readr)
library(dplyr)
library(data.table)
library(ggplot2)
library(doBy)
library(reshape2)

# Read in sample snp count files (script: calculate_pn_ps_part2.R)
snps <- read_delim("pop_level_number_syn_nonsyn_SNPs.txt", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE)

# Read in total syn and non syn possible given the genome (script: calculate_pn_ps_parti.R)
avaliable_syn_nonsyn <- read_delim("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/avaliable_syn_nonsyn.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

# Merge the above to each sample and calculate the proportion of syn and nonsyn
# Add in the actual ratio then of pn/ps
# Title the column with the sample name
snps <- merge(snps, avaliable_syn_nonsyn, by = "gene_id")
snps$pn <- snps$non_syn/snps$non_syn_all_gene
snps$ps <- snps$syn/snps$syn_all_gene
snps$pn_ps <- snps$pn/snps$ps
snps <- snps[,c(2,1,9)]

snps[sapply(snps, is.infinite)] <- NA

write.table(snps, file="pn_ps_population_level.txt", sep="\t", quote = F, col.names = T, row.names = F)
all_data <- snps

# Have a look at some summary stats
head(all_data)
hist(all_data$pn_ps)
nrow(all_data) #8563
all_data <- all_data[!is.na(all_data$pn_ps),] #7678
nrow(all_data[all_data$pn_ps >1,]) #199

# Make a plot
head(all_data)
all_data$chr[all_data$chr %like% "NW"] <- "unplaced_scaffold"
ggplot(all_data, aes(y=pn_ps, x=chr))+
  geom_boxplot()

# Correlate it with methylation
setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/Previous/Differential_methylation/weighted_meth/just_male_queen")
file.list = list.files(("./"),pattern="*meth_all_features.txt")

read_filei <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_filei)
sample_names <- list("m08","m19","m23","m37","q08","q19","q23","q37")
names(samples) <- sample_names

for(i in seq_along(samples)){
  samples[[i]] <- samples[[i]][samples[[i]]$feature=="gene",]
  samples[[i]] <- samples[[i]][,c(2,9)]
  samples[[i]]$sample <- names(samples[i])
}

meth_data <- bind_rows(samples)
meth_data$sex <- "male"
meth_data$sex[meth_data$sample %like% "q"] <- "queen"
table(meth_data$sex)
meth_data <- meth_data[,-3]
meth_data2 <- dcast(meth_data, gene_id ~ sex, value.var="weightedMeth", fun.aggregate = mean)

both <- merge(all_data, meth_data2, by = c("gene_id")) #7676
both <- both[complete.cases(both),] #7666
both$weightedMeth_mean <- (both$male + both$queen)/2

head(both)
ggplot(both, aes(x=weightedMeth_mean, y=pn_ps))+
  geom_point()

both$selection <- "< 1"
both$selection[both$pn_ps>1] <- "> 1"

both$logmeth <- log10(both$weightedMeth_mean)

ggplot(both, aes(x=selection, y=logmeth))+
  geom_violin()+
  #geom_boxplot()+
  #geom_point(position = "jitter")+
  xlab("pN / pS")+
  ylab("log10(Weighted Methylation)")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())

# stats?
head(both)

both$meth_category <- "methylated"
both$meth_category[both$weightedMeth_mean < 0.005] <- "unmethylated"
table(both$meth_category)


ggplot(both, aes(x=meth_category, y=pn_ps))+
  geom_boxplot()+
  xlab("Gene Category")+
  ylab("pN / pS")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())
  


# break down by sex and high, med, low, no meth status
head(both)

by_sex <- both[,c(1,2,3,4,5,7)]
head(by_sex)
by_sex_long <- melt(by_sex, id.vars = c("gene_id","chr","pn_ps","selection"))
colnames(by_sex_long) <- c("gene_id","chr","pn_ps","selection","sex","methylation")
by_sex_long$logmeth <- log10(by_sex_long$methylation)

ggplot(by_sex_long, aes(x=selection, y=logmeth, fill=sex))+
  geom_violin()+
  #geom_boxplot()+
  #geom_point(position = "jitter")+
  xlab("pN / pS")+
  ylab("log10(Weighted Methylation)")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("male","queen"),
                    values = c("midnightblue","#CC6677"),
                    labels = c("Male","Queen"))

library(multcomp)
head(by_sex_long)
table(by_sex_long$selection)
model1<-lm(pn_ps ~ methylation * sex , data=by_sex_long)
model2<-lm(pn_ps ~ methylation  , data=by_sex_long)
anova(model1,model2) # No sig interaction
summary.lm(model2) # Not sig



# also diff meth genes

# Next question what if we do this by highly methylated genes and lowly methylated genes?
# Also is there a relationship with differentially methylated genes?