# Edit coverage files to calculate meth per gene levels

gunzip *cov.gz

for file in $(ls *cov)
do
    base=$(basename ${file} "_1_bismark_bt2_pe.deduplicated.bismark.cov")
    cut -f1,2,5,6 ${file} > ${base}_coverage.txt
done

# --- R ---
module load R/3.6.1

library(readr)

file.list = list.files(("./"),pattern="*_coverage.txt")

read_file1 <- function(x){
  read_delim(x, "\t", col_names=F)
}

samples <- lapply(file.list, read_file1)

sample_names <- list("m08","m19","m23","m37","q08","q19","q23","q37","w08","w19","w23","w37")
names(samples) <- sample_names

for(i in seq_along(samples)){
    colnames(samples[[i]]) <- c("chr", "cpg", "count_c", "count_t")
    samples[[i]]$total_coverage <- samples[[i]]$count_c + samples[[i]]$count_t
    samples[[i]] <- samples[[i]][,-4]
    final_file <- samples[[i]]
    myfile <- file.path("./", paste0(names(samples[i]),"_","final_coverage.txt"))
    write.table(final_file, file=myfile, quote=F, sep="\t", row.names=F)
}