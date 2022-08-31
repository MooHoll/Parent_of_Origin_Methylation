#-----------------------------------------------
# Script created by Alun Jones, see paper Bebane et al. (2019) Neonics and bumblebees...
#-----------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/PoO_Methylation_BB/New_2022/GO_analysis")
library(GOstats)
library(GSEABase)
library(treemap)

#-----------------------------------------------
# Read in background GO set and make compatible with GOstats

#GO_annotations <- read.table("./GO_male_all_methylated_genes.txt")
#GO_annotations <- read.table("./GO_queen_all_methylated_genes.txt")
#GO_annotations <- read.table("./GO_worker_all_methylated_genes.txt")

#GO_annotations <- read.table("./MvW_all_background.txt")
#GO_annotations <- read.table("./MvQ_all_background.txt")
#GO_annotations <- read.table("./WvQ_all_background.txt")

#GO_annotations <- read.table("./GO_diff_meth_MvW.txt")
#GO_annotations <- read.table("./GO_diff_meth_MvQ.txt")
#GO_annotations <- read.table("./GO_diff_meth_WvQ.txt")

GO_annotations <- read.table("./all_comparisons_background.txt")

#-----------------------------------------------
GO_annotations[,3] <- paste("IEA")
names(GO_annotations) <- c("genes","GOIds","evi")
GO_annotations[,3] <- paste("IEA")
GO_annotations <- GO_annotations[c(2,3,1)]
GO_annotations <- GO_annotations[-1,] # remove the original colnames

#-----------------------------------------------
# Create necessary objects

GO_frame <- GOFrame(GO_annotations,organism = "Bombus terrestris")
goAllFrame <- GOAllFrame(GO_frame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
universe <- as.vector(unique(GO_annotations[,3]))

#-----------------------------------------------
# Read in gene's of interest 

#my_genes <- read.table("./male_high_methylated_genes.txt", header = T)
#my_genes <- read.table("./queen_high_methylated_genes.txt", header = T)
#my_genes <- read.table("./worker_high_methylated_genes.txt", header = T)

#my_genes <- read.table("./MvW_all_diff_meth_genes.txt", header = T)
#my_genes <- read.table("./MvQ_all_diff_meth_genes.txt", header = T)
#my_genes <- read.table("./WvQ_all_diff_meth_genes.txt", header = T)

#my_genes <- read.table("./MvW_male_hypermeth_genes.txt", header = T)
#my_genes <- read.table("./MvW_worker_hypermeth_genes.txt", header = T)

#my_genes <- read.table("./MvQ_male_hypermeth_genes.txt", header = T)
#my_genes <- read.table("./MvQ_queen_hypermeth_genes.txt", header = T)

#my_genes <- read.table("./WvQ_worker_hypermeth_genes.txt", header = T)
#my_genes <- read.table("./WvQ_queen_hypermeth_genes.txt", header = T)

#my_genes <- read.table("./GO_diff_meth_hyper_WQvM.txt", header = T)
my_genes <- read.table("./GO_diff_meth_hyper_MvWQ.txt", header = T)

#-----------------------------------------------
my_genes <- as.data.frame(na.omit(my_genes$gene_id))
colnames(my_genes) <- "genes"
my_genes <- as.vector(my_genes[,1])

#-----------------------------------------------
# Keep only genes with annotated GOs

my_genes <- my_genes[my_genes %in% universe]
length(my_genes)

# High meth vs all meth for indiv sexes and castes
# Male: 27/37
# Queen: 12/18
# Worker 15/21

# Diff meth vs all meth found in both in the comparison
# MvW: 150/161
# MvQ: 147/161
# WvQ: 55/59

# Hyper meth vs diff meth per comparison
# Male (in MvW): 91/97
# Worker (in MvW): 93/99
# Male (in MvQ): 91/96
# Queen (in MvQ): 102/111
# Worker (in WvQ): 38/39
# Queen (in WvQ): 29/33

# QW common vs male with all diff meth as background
# WQ hyper: 24/28
# M hyper: 22/23

#-----------------------------------------------
# Set up paramters for hypergeometric test

Get_GO_params_all <- function(genes_of_i,universe,pvalue_cut){
  onto_terms <- c("BP","CC","MF")
  directions <- c("over","under")
  param_list <- list()
  name_1 <- list()
  for(i in 1:3){
    for(j in 1:2){
      name_1 <- c(name_1,paste(onto_terms[i],directions[j],sep = "_"))
      parameters <- GSEAGOHyperGParams(name="Hygeo params",
                                       geneSetCollection = gsc,
                                       universeGeneIds = universe,
                                       geneIds = my_genes,
                                       ontology = paste(onto_terms[i]),
                                       pvalueCutoff = pvalue_cut,
                                       conditional = T,testDirection = paste(directions[j]))
      param_list <- c(param_list,parameters)
    }
  }
  names(param_list) <- name_1
  return(param_list)
}

param_list <- Get_GO_params_all(genes_of_i = DE_Genes_A,universe = universe,
                                pvalue_cut = 0.05)


#-----------------------------------------------
# Hypergeometric test

Hyper_G_test <- function(param_list){
  Hyper_G_list <- list()
  for(i in 1:length(param_list)){
    res <- hyperGTest(param_list[[i]])
    Hyper_G_list <- c(Hyper_G_list,res)
  }
  names(Hyper_G_list) <- names(param_list)
  return(Hyper_G_list)
}


GO_enrichment <- Hyper_G_test(param_list = param_list)


#-----------------------------------------------
# Choose what you want to look for, here the choice is biological process over-represented 

Result <- summary(GO_enrichment[["BP_over"]])

#-----------------------------------------------
# FDR correction

Result_FDR <- Result[p.adjust(Result$Pvalue,method = "fdr") < 0.05,]

#-----------------------------------------------
# Make an output for REVIGO and write out

REVIGO <- Result_FDR[,1:2]

#write.table(REVIGO,"./outputs/male_high_meth_against_all_meth.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./outputs/queen_high_meth_against_all_meth.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./outputs/worker_high_meth_against_all_meth.txt",row.names = F,sep = "\t",quote = F)

#write.table(REVIGO,"./outputs/all_diff_meth_MvW.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./outputs/all_diff_meth_MvQ.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./outputs/all_diff_meth_WvQ.txt",row.names = F,sep = "\t",quote = F)

#write.table(REVIGO,"./outputs/male_hyper_MvW.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./outputs/worker_hyper_MvW.txt",row.names = F,sep = "\t",quote = F)

#write.table(REVIGO,"./outputs/male_hyper_MvQ.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./outputs/queen_hyper_MvQ.txt",row.names = F,sep = "\t",quote = F)

#write.table(REVIGO,"./outputs/worker_hyper_WvQ.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./outputs/queen_hyper_WvQ.txt",row.names = F,sep = "\t",quote = F)

#write.table(REVIGO,"./outputs/WQ_hyper_vsM.txt",row.names = F,sep = "\t",quote = F)
write.table(REVIGO,"./outputs/M_hyper_vsWQ.txt",row.names = F,sep = "\t",quote = F)



