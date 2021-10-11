# Making the file which has the total CpGs per gene information

library(sqldf)
library(readr)
library(doBy)

cpgs <- read.delim("destranded_cpg_positions.txt", header=F)
genes <- read.delim("annotationfile_with_promotors.csv", header=T)

output <- sqldf("SELECT sg.V1,
                sg.V2,
                fg.chr,
                fg.start,
                fg.end,
                fg.feature,
                fg.geneID
                FROM cpgs AS sg
                LEFT JOIN genes AS fg 
                ON sg.V1 = fg.chr
                AND (sg.V2 >= fg.start AND sg.V2 <= fg.end)")
output <- output[!is.na(output$geneID),]
nrow(output) 

output$cpg_counter <- 1

final <- summaryBy(cpg_counter ~ geneID+chr+start+end, data = output, FUN=sum)
nrow (final) 

write.table(final, file="annotation_with_total_cpgs.txt", col.names=T, row.names=F, quote=F)
