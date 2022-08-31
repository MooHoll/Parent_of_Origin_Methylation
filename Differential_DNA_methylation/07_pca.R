## -------------------------------------------------------------------------
## Making Fancy Genome-Methylation Differences Graphs
## -------------------------------------------------------------------------

# Load packages etc.
setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/coverage_files")
library(grid)
library(readr)
library(ggplot2)
library(methylKit)
library(ggdendro)
library(ggfortify)
library(ggrepel)

## -------------------------------------------------------------------------

# Get a methylkit object for all samples
sample.list <- list("m08_merged_CpG_evidence.cov" ,"m19_merged_CpG_evidence.cov",
                    "m23_merged_CpG_evidence.cov","m37_merged_CpG_evidence.cov",
                    "q08_merged_CpG_evidence.cov" ,"q19_merged_CpG_evidence.cov",
                    "q23_merged_CpG_evidence.cov","q37_merged_CpG_evidence.cov",
                    "w08_merged_CpG_evidence.cov" ,"w19_merged_CpG_evidence.cov",
                    "w23_merged_CpG_evidence.cov","w37_merged_CpG_evidence.cov")

CPGRaw <- methRead(sample.list, 
                   sample.id = list("m08", "m19","m23","m37",
                                    "q08", "q19","q23","q37",
                                    "w08", "w19","w23","w37"),
                   assembly="bter_1.0",
                   treatment=c(0,0,0,0,1,1,1,1,2,2,2,2),
                   context="CpG",
                   dbtype = NA,
                   pipeline = "bismarkCoverage",
                   header = T, sep = "\t", mincov=1,
                   dbdir= path)

filtered_data <- filterByCoverage(CPGRaw,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)
rm(CPGRaw)
normalized <- normalizeCoverage(filtered_data)
rm(filtered_data)
meth_all_data <- unite(normalized, destrand=F) 
rm(normalized)
nrow(meth_all_data) # 1463570

df_meth_all <- getData(meth_all_data)

a <- df_meth_all[,5:6]
b <- df_meth_all[,8:9]
c <- df_meth_all[,11:12]
d <- df_meth_all[,14:15]
e <- df_meth_all[,17:18]
f <- df_meth_all[,20:21]
g <- df_meth_all[,23:24]
h <- df_meth_all[,26:27]
k <- df_meth_all[,29:30]
l <- df_meth_all[,32:33]
m <- df_meth_all[,35:36]
n <- df_meth_all[,38:39]

bt <- function(a, b, p = 0.005) {binom.test(a, b, 0.005, alternative="greater")$p.value}
allrows <- data.frame(CT=numeric(),Ccount=numeric(),pVal=numeric(),FDR=numeric(),row=numeric())

for (df in list(a,b,c,d,e,f,g,h,k,l,m,n)) {
  colnames(df) <- c("CT", "Ccount")
  df$pVal <- mapply(bt, df$Ccount, df$CT)
  df$FDR <- p.adjust(df$pVal, method = "BH", n = length(df$pVal))
  df$row <- as.numeric(rownames(df))
  dfmeth <- subset(df, df$FDR < 0.05)
  allrows <- rbind(allrows, dfmeth)
}

meth_positions <- as.vector(as.numeric(unique(allrows$row))) 
length(meth_positions) # 6514

subset_methBase <- methylKit::select(meth_all_data, meth_positions)
head(subset_methBase)

## -------------------------------------------------------------------------
# Make a nice PCA

PCA_data <- PCASamples(subset_methBase, screeplot=F, obj.return = T)
PCA_data1 <- as.data.frame(PCA_data$x)
PCA_data1$sample <- row.names(PCA_data1)

PCA_data1$Caste <- c(rep("Male", 4), rep("Queen", 4), rep("Worker", 4))


percentage <- round(PCA_data$sdev / sum(PCA_data$sdev) * 100, 0)
percentage <- paste(colnames(PCA_data), "(", paste( as.character(percentage), "%", ")", sep="") )


ggplot(PCA_data1, aes(PC1, PC2, colour=Caste))+
  geom_point(size=14)+
  geom_text_repel(aes(label=sample), size=12,show.legend=FALSE, 
                  point.padding = 2, box.padding = 1)+
  theme_bw()+
  xlab(paste0("PC1:",percentage[1],"variance")) +
  ylab(paste0("PC2:",percentage[2],"variance")) +
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=30),
        legend.text=element_text(size=30),
        legend.title=element_blank())+
  scale_colour_manual(values=c("#44AA99","#CC6677","#DDCC77"))
