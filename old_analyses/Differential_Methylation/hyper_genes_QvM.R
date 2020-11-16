# Have a look at hypermethylated genes in either male/queen

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Methylation/Diff_Meth/queens_vs_males")

library(readr)

all_data <- read_csv("DMGs_MSCfilter_with_geneID_Male_Queen.csv")
check <- as.data.frame(unique(all_data$geneID)) #412
#write_delim(check, "diffmeth_unique_genes.txt")

# Methylkit works as 1-0 and males =1 and queens =0, so +ve value is hyper in male
male_hyper <- subset(all_data, meth_diff>0) 
male_hyper <- as.data.frame(unique(male_hyper$geneID)) #211 - 49 in common = 162
colnames(male_hyper) <-"geneID"
#write_delim(male_hyper, "male_hyper_genes.txt")

queen_hyper <- subset(all_data, meth_diff<0) 
queen_hyper <- as.data.frame(unique(queen_hyper$geneID)) #250 - 49 = 201
colnames(queen_hyper) <-"geneID"
#write_delim(queen_hyper, "queen_hyper_genes.txt")

#This is a goodness of fit
# Testing if more hyper meth in queen or male
observed = c(211, 250)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions

chisq.test(x = observed,
           p = expected)
# X-squared = 3.2993, df = 1, p-value = 0.06931


# Check if any genes contains a CpG hyper in male and queen
both <- merge(male_hyper, queen_hyper, by="geneID") #49
#write_delim(both, "both_hyper_genes.txt")


## Working out weighted meth so can do the 10% difference threshold over the whole gene

library(reshape2)

weighted_meth_by_sample_type <- read_delim("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Methylation/PoO_meth/scatter_graphs/weighted_meth/weighted_meth_by_sample_type.txt", 
                                           "\t", escape_double = FALSE, trim_ws = TRUE)
queen_info <- subset(weighted_meth_by_sample_type, origin == "queens")
colnames(queen_info)[3] <- "queens"
queen_info <- queen_info[,-2]

male_info <- subset(weighted_meth_by_sample_type, origin == "male")
colnames(male_info)[3] <- "males"
male_info <- male_info[,-2]

both_data <- merge(male_info, queen_info, by = "geneID")
both_data$diff <- 100*((both_data$males - both_data$queens)/both_data$queens)

both_data <- both_data[!is.na(both_data$diff),]
both_data$diff[both_data$diff==Inf] <- 100
filtered_data <- both_data[(both_data$diff > 10 | both_data$diff < -10),]

#344 with 10% weighted meth diff
all_data_look <- filtered_data[(filtered_data$geneID %in% check$`unique(all_data$geneID)`),]
all_data_look$hypermethylated <- "queen"
all_data_look$hypermethylated[all_data_look$diff > 0] <- "male"
write.table(all_data_look, file="diff_meth_QvM_10%_diff_all_info.txt", sep="\t", quote = F, col.names = T, row.names = F )


genes <- as.data.frame(all_data_look$geneID)
colnames(genes) <- "geneID"
write.table(genes, file="diff_meth_weighted_genes.txt", sep="\t", quote = F, col.names = T, row.names = F )



male_hyper <- subset(all_data_look, diff > 0) # 151
queen_hyper <- subset(all_data_look, diff < 0) # 193

# Goodness of fit
observed = c(151, 193)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions

chisq.test(x = observed,
           p = expected) # X-squared = 5.1279, df = 1, p-value = 0.02354

male_hyper <- as.data.frame(male_hyper$geneID)
colnames(male_hyper) <- "geneID"
queen_hyper <- as.data.frame(queen_hyper$geneID)
colnames(queen_hyper) <- "geneID"

write.table(male_hyper, file="male_hyper_genes.txt", sep="\t", quote = F, col.names = T, row.names = F )
write.table(queen_hyper, file="queen_hyper_genes.txt", sep="\t", quote = F, col.names = T, row.names = F )


