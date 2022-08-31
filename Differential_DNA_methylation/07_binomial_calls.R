## -------------------------------------------------------------------------
## Methylation binomial calls
## -------------------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/cov_reports")
library(methylKit)
library(readr)
library(dplyr)

# -------------------------------------------------------------------------

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

a <- df_meth_all[,c(1:2,5:6)]
b <- df_meth_all[,c(1:2,8:9)]
c <- df_meth_all[,c(1:2,11:12)]
d <- df_meth_all[,c(1:2,14:15)]
e <- df_meth_all[,c(1:2,17:18)]
f <- df_meth_all[,c(1:2,20:21)]
g <- df_meth_all[,c(1:2,23:24)]
h <- df_meth_all[,c(1:2,26:27)]
k <- df_meth_all[,c(1:2,29:30)]
l <- df_meth_all[,c(1:2,32:33)]
m <- df_meth_all[,c(1:2,35:36)]
n <- df_meth_all[,c(1:2,38:39)]

datalist = list(a,b,c,d,e,f,g,h,k,l,m,n)

sample.id = list("m08", "m19","m23","m37",
                 "q08", "q19","q23","q37",
                 "w08", "w19","w23","w37")
names(datalist) <- sample.id

# Probability of success = average non-conversion rate
bt <- function(a, b, p = 0.0036) {binom.test(a, b, 0.0036, alternative="greater") $p.value}

for (i in seq_along(datalist)) {  
  colnames(datalist[[i]]) <- c("chr","position","coverage", "Ccount")
  datalist[[i]] <- datalist[[i]][!is.na(datalist[[i]]$coverage),]
  datalist[[i]]$pVal <- mapply(bt, datalist[[i]]$Ccount, datalist[[i]]$coverage)
  datalist[[i]]$FDR <- p.adjust(datalist[[i]]$pVal, method = "BH", n = length(datalist[[i]]$pVal))
  datalist[[i]]$methylated <- 0
  datalist[[i]]$methylated[datalist[[i]]$FDR < 0.05] <- 1
  datalist[[i]] <- datalist[[i]][,c(1,2,7)]
  colnames(datalist[[i]]) <- c("chr","position",paste0(names(datalist[i])))
}

all = Reduce(function(...) merge(..., all=T, by = c("chr","position")), datalist)
write.table(all, file="methylation_calls_per_sample.txt", quote=F, sep="\t", row.names=F, col.names=T)

