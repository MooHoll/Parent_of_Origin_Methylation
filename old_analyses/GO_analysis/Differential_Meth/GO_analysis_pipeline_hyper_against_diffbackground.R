
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Methylation/Diff_Meth/queens_vs_males/GO_analysis")
################################################
### GO ANALYSIS packages

### loading up the packages

if (!require("GOstats")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("GOstats")
  library(GOstats)
}

if (!require("GSEABase")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("GSEABase")
  library(GSEABase)
}

if (!require("treemap")) {
  install.packages("treemap", dependencies = TRUE)
  library(treemap)
}

################################################
### reading in the GO terms
### and sorting out for GOStats

GO_annotations <- read.table("GOs_all_diffmeth_genes_as_background.txt")

GO_annotations[,3] <- paste("IEA")

names(GO_annotations) <- c("genes","GOIds","evi")

GO_annotations[,3] <- paste("IEA")

GO_annotations <- GO_annotations[c(2,3,1)]


################################################
### creating a GO frame object and the 
### gene set collection object
### creating gene universe

GO_frame <- GOFrame(GO_annotations,organism = "Bombus terrestris")

goAllFrame <- GOAllFrame(GO_frame)

gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

universe <- as.vector(unique(GO_annotations[,3]))


################################################
### getting the list of genes of interest eg
### differentially expressed Loci

# Change file here depending on gene list
#my_genes <- read.csv("./gene_lists/queen_hyper_genes.txt",
 #                    header=T)
my_genes <- read.csv("./gene_lists/male_hyper_genes.txt",
                     header=T)
my_genes <- as.data.frame(na.omit(my_genes$geneID))
colnames(my_genes) <- "genes"

# change class to a vector
my_genes <- as.vector(my_genes[,1])


# Keeping only DE genes that been annotated
my_genes <- my_genes[my_genes %in% universe]
#174/193 (hyper_queen) 
#134/151 (hyper_male) 

################################################
### Setting up the parameters for the 
### hypergeometric test
### function which produces list of different types of go analysis 
### eg Molecular function over and under respresentes GO terms

Get_GO_params_all <- function(genes_of_i,universe,pvalue_cut){
  
  onto_terms <- c("BP","CC","MF")
  directions <- c("over","under")
  param_list <- list()
  name_1 <- list()
  for(i in 1:3){
    for(j in 1:2){
      name_1 <- c(name_1,paste(onto_terms[i],directions[j],sep = "_"))
      parameters <- GSEAGOHyperGParams(name="Bumblebee Hygeo params",
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





################################################
### conducting the hypergeometric test
### using the paramters set
### produces a list of results from the list of paramters


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

### P adjustment value to incorporate with p.adjust function

################################################
### in order to get the desired GO enrichment
### use summary function on element of
### list which is the enrichment you want to
### look at eg Molecular function overrepresentation
### also using a fdr adjusment.  Can use other adjusments
### but because they so stringent reliable to have no significant terms


# Here the choice is biological process over-represented (checked online and this is best choice)
Result <- summary(GO_enrichment[["BP_over"]])

Result_FDR <- Result[p.adjust(Result$Pvalue,method = "fdr") < 0.05,]

### take first two coloumns into revigo website

REVIGO <- Result_FDR[,1:2]

#write.table(REVIGO,"enriched_GOs_hyper_queen_against_diffmeth.txt",row.names = F,sep = "\t",quote = F)
write.table(REVIGO,"enriched_GOs_hyper_male_against_diffmeth.txt",row.names = F,sep = "\t",quote = F)
