## -------------------------------------------------------------------------
## Weighted Methylation per annotation: make file with cpg count per annoation
## -------------------------------------------------------------------------

# Load packages etc.
library(sqldf)
library(readr)
library(doBy)
library(foreach)
library(doParallel)

registerDoParallel(cores=4)

## -------------------------------------------------------------------------
# Making the file which has the total CpGs per gene information

file.list = list.files(("./"),pattern="*_total_cpgs_in_genome.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = F, trim_ws = T)
}

cpgs <- lapply(file.list, read_file1)
sample_names <- list("m08","m19","m23","m37")
names(cpgs) <- sample_names

for(i in seq_along(cpgs)){
  colnames(cpgs[[i]]) <- c("chr", "cpg_position")
  cpgs[[i]]$cpg_position <- as.numeric(cpgs[[i]]$cpg_position)
}


# --------------------------------------------------------------------

genes <- read.delim("ref_Bter_1.0_top_level_numbered_exons.txt", header=T)
genes$start <- as.numeric(genes$start)
genes$end <- as.numeric(genes$end)

foreach(i = seq_along(cpgs)) %dopar% {
  df <- cpgs[[i]]
  output <- sqldf("SELECT cpgs.chr,
                cpgs.cpg_position,
                genes.chr,
                genes.feature,
                genes.start,
                genes.end,
                genes.gene_id,
                genes.number
                FROM df AS cpgs
                LEFT JOIN genes AS genes 
                ON cpgs.chr = genes.chr
                AND (cpgs.cpg_position >= genes.start AND cpgs.cpg_position <= genes.end)")
  output <- output[!is.na(output$gene_id),]
  output$cpg_counter <- 1
  final <- summaryBy(cpg_counter ~ gene_id+chr+start+end+feature+number, data = output, FUN=sum)
  myfile <- file.path("./", paste0(names(cpgs[i]),"_","annotation_with_total_cpgs.txt"))
  write.table(final, file=myfile, col.names=T, row.names=F, quote=F)
}