### Overlapping geneID with SNP positions using a databse 

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Count_Table")

library(dplyr)
library(sqldf)

# Main file from GFF containing gene with start and end position
final_genes_sql <- read.csv.sql("./all_genes_from_genome.txt",
                       sql = "select * from file", sep = "\t", header = FALSE) 

# Testing with just one file 
sample_genes_sql <- read.csv.sql("./read_counts/M_trimmed_02_10_KOP_male_counts_only.txt",
                                sql = "select * from file", sep = "\t", header = FALSE) 

output <- sqldf("SELECT sg.V1,
      sg.V2,
      sg.V3,
      fg.V4
      FROM sample_genes_sql AS sg
      LEFT JOIN final_genes_sql AS fg 
      ON sg.V1 = fg.V1
      AND (sg.V2 >= fg.V2 AND sg.V2 <= fg.V3)")

# Removing rows that didn't match, checked output with Vagelis's Python script and takes 2mins per file
# compared to 11hours per file
m0210_male<-subset(output, !output$V4=="NA")

### Need to rbind all samples for male and female SNPs, so 256 files goes to 128
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Count_Table/read_counts")

file.list = list.files(pattern="*.txt")

for (i in file.list){
  x<-read.csv.sql(i, header=F, sep="\t")
  assign(i,x)
}

m0210abd<-rbind(M_trimmed_02_10_ABD_male_counts_only.txt, M_trimmed_02_10_ABD_queen_counts_only.txt)
m0212abd<-rbind(M_trimmed_02_12_ABD_male_counts_only.txt, M_trimmed_02_12_ABD_queen_counts_only.txt)
m0218abd<-rbind(M_trimmed_02_18_ABD_male_counts_only.txt, M_trimmed_02_18_ABD_queen_counts_only.txt)
m0222abd<-rbind(M_trimmed_02_22_ABD_male_counts_only.txt, M_trimmed_02_22_ABD_queen_counts_only.txt)
m0235abd<-rbind(M_trimmed_02_35_ABD_male_counts_only.txt, M_trimmed_02_35_ABD_queen_counts_only.txt)
m0237abd<-rbind(M_trimmed_02_37_ABD_male_counts_only.txt, M_trimmed_02_37_ABD_queen_counts_only.txt)
m0244abd<-rbind(M_trimmed_02_44_ABD_male_counts_only.txt, M_trimmed_02_44_ABD_queen_counts_only.txt)
m0254abd<-rbind(M_trimmed_02_54_ABD_male_counts_only.txt, M_trimmed_02_54_ABD_queen_counts_only.txt)
m1274abd<-rbind(M_trimmed_12_74_ABD_male_counts_only.txt, M_trimmed_12_74_ABD_queen_counts_only.txt)
m1275abd<-rbind(M_trimmed_12_75_ABD_male_counts_only.txt, M_trimmed_12_75_ABD_queen_counts_only.txt)
m1276abd<-rbind(M_trimmed_12_76_ABD_male_counts_only.txt, M_trimmed_12_76_ABD_queen_counts_only.txt)
m1278abd<-rbind(M_trimmed_12_78_ABD_male_counts_only.txt, M_trimmed_12_78_ABD_queen_counts_only.txt)
m1282abd<-rbind(M_trimmed_12_82_ABD_male_counts_only.txt, M_trimmed_12_82_ABD_queen_counts_only.txt)
m1284abd<-rbind(M_trimmed_12_84_ABD_male_counts_only.txt, M_trimmed_12_84_ABD_queen_counts_only.txt)
m1285abd<-rbind(M_trimmed_12_85_ABD_male_counts_only.txt, M_trimmed_12_85_ABD_queen_counts_only.txt)
m1289abd<-rbind(M_trimmed_12_89_ABD_male_counts_only.txt, M_trimmed_12_89_ABD_queen_counts_only.txt)
m2209abd<-rbind(M_trimmed_22_09_ABD_male_counts_only.txt, M_trimmed_22_09_ABD_queen_counts_only.txt)
m2218abd<-rbind(M_trimmed_22_18_ABD_male_counts_only.txt, M_trimmed_22_18_ABD_queen_counts_only.txt)
m2221abd<-rbind(M_trimmed_22_21_ABD_male_counts_only.txt, M_trimmed_22_21_ABD_queen_counts_only.txt)
m2223abd<-rbind(M_trimmed_22_23_ABD_male_counts_only.txt, M_trimmed_22_23_ABD_queen_counts_only.txt)
m2226abd<-rbind(M_trimmed_22_26_ABD_male_counts_only.txt, M_trimmed_22_26_ABD_queen_counts_only.txt)
m2230abd<-rbind(M_trimmed_22_30_ABD_male_counts_only.txt, M_trimmed_22_30_ABD_queen_counts_only.txt)
m2235abd<-rbind(M_trimmed_22_35_ABD_male_counts_only.txt, M_trimmed_22_35_ABD_queen_counts_only.txt)
m2249abd<-rbind(M_trimmed_22_49_ABD_male_counts_only.txt, M_trimmed_22_49_ABD_queen_counts_only.txt)
m3144abd<-rbind(M_trimmed_31_44_ABD_male_counts_only.txt, M_trimmed_31_44_ABD_queen_counts_only.txt)
m3145abd<-rbind(M_trimmed_31_45_ABD_male_counts_only.txt, M_trimmed_31_45_ABD_queen_counts_only.txt)
m3164abd<-rbind(M_trimmed_31_64_ABD_male_counts_only.txt, M_trimmed_31_64_ABD_queen_counts_only.txt)
m3165abd<-rbind(M_trimmed_31_65_ABD_male_counts_only.txt, M_trimmed_31_65_ABD_queen_counts_only.txt)
m3179abd<-rbind(M_trimmed_31_79_ABD_male_counts_only.txt, M_trimmed_31_79_ABD_queen_counts_only.txt)
m3184abd<-rbind(M_trimmed_31_84_ABD_male_counts_only.txt, M_trimmed_31_84_ABD_queen_counts_only.txt)
m3186abd<-rbind(M_trimmed_31_86_ABD_male_counts_only.txt, M_trimmed_31_86_ABD_queen_counts_only.txt)
m3189abd<-rbind(M_trimmed_31_89_ABD_male_counts_only.txt, M_trimmed_31_89_ABD_queen_counts_only.txt)

m0210kop<-rbind(M_trimmed_02_10_KOP_male_counts_only.txt, M_trimmed_02_10_KOP_queen_counts_only.txt)
m0212kop<-rbind(M_trimmed_02_12_KOP_male_counts_only.txt, M_trimmed_02_12_KOP_queen_counts_only.txt)
m0218kop<-rbind(M_trimmed_02_18_KOP_male_counts_only.txt, M_trimmed_02_18_KOP_queen_counts_only.txt)
m0222kop<-rbind(M_trimmed_02_22_KOP_male_counts_only.txt, M_trimmed_02_22_KOP_queen_counts_only.txt)
m0235kop<-rbind(M_trimmed_02_35_KOP_male_counts_only.txt, M_trimmed_02_35_KOP_queen_counts_only.txt)
m0237kop<-rbind(M_trimmed_02_37_KOP_male_counts_only.txt, M_trimmed_02_37_KOP_queen_counts_only.txt)
m0244kop<-rbind(M_trimmed_02_44_KOP_male_counts_only.txt, M_trimmed_02_44_KOP_queen_counts_only.txt)
m0254kop<-rbind(M_trimmed_02_54_KOP_male_counts_only.txt, M_trimmed_02_54_KOP_queen_counts_only.txt)
m1274kop<-rbind(M_trimmed_12_74_KOP_male_counts_only.txt, M_trimmed_12_74_KOP_queen_counts_only.txt)
m1275kop<-rbind(M_trimmed_12_75_KOP_male_counts_only.txt, M_trimmed_12_75_KOP_queen_counts_only.txt)
m1276kop<-rbind(M_trimmed_12_76_KOP_male_counts_only.txt, M_trimmed_12_76_KOP_queen_counts_only.txt)
m1278kop<-rbind(M_trimmed_12_78_KOP_male_counts_only.txt, M_trimmed_12_78_KOP_queen_counts_only.txt)
m1282kop<-rbind(M_trimmed_12_82_KOP_male_counts_only.txt, M_trimmed_12_82_KOP_queen_counts_only.txt)
m1284kop<-rbind(M_trimmed_12_84_KOP_male_counts_only.txt, M_trimmed_12_84_KOP_queen_counts_only.txt)
m1285kop<-rbind(M_trimmed_12_85_KOP_male_counts_only.txt, M_trimmed_12_85_KOP_queen_counts_only.txt)
m1289kop<-rbind(M_trimmed_12_89_KOP_male_counts_only.txt, M_trimmed_12_89_KOP_queen_counts_only.txt)
m2209kop<-rbind(M_trimmed_22_09_KOP_male_counts_only.txt, M_trimmed_22_09_KOP_queen_counts_only.txt)
m2218kop<-rbind(M_trimmed_22_18_KOP_male_counts_only.txt, M_trimmed_22_18_KOP_queen_counts_only.txt)
m2221kop<-rbind(M_trimmed_22_21_KOP_male_counts_only.txt, M_trimmed_22_21_KOP_queen_counts_only.txt)
m2223kop<-rbind(M_trimmed_22_23_KOP_male_counts_only.txt, M_trimmed_22_23_KOP_queen_counts_only.txt)
m2226kop<-rbind(M_trimmed_22_26_KOP_male_counts_only.txt, M_trimmed_22_26_KOP_queen_counts_only.txt)
m2230kop<-rbind(M_trimmed_22_30_KOP_male_counts_only.txt, M_trimmed_22_30_KOP_queen_counts_only.txt)
m2235kop<-rbind(M_trimmed_22_35_KOP_male_counts_only.txt, M_trimmed_22_35_KOP_queen_counts_only.txt)
m2249kop<-rbind(M_trimmed_22_49_KOP_male_counts_only.txt, M_trimmed_22_49_KOP_queen_counts_only.txt)
m3144kop<-rbind(M_trimmed_31_44_KOP_male_counts_only.txt, M_trimmed_31_44_KOP_queen_counts_only.txt)
m3145kop<-rbind(M_trimmed_31_45_KOP_male_counts_only.txt, M_trimmed_31_45_KOP_queen_counts_only.txt)
m3164kop<-rbind(M_trimmed_31_64_KOP_male_counts_only.txt, M_trimmed_31_64_KOP_queen_counts_only.txt)
m3165kop<-rbind(M_trimmed_31_65_KOP_male_counts_only.txt, M_trimmed_31_65_KOP_queen_counts_only.txt)
m3179kop<-rbind(M_trimmed_31_79_KOP_male_counts_only.txt, M_trimmed_31_79_KOP_queen_counts_only.txt)
m3184kop<-rbind(M_trimmed_31_84_KOP_male_counts_only.txt, M_trimmed_31_84_KOP_queen_counts_only.txt)
m3186kop<-rbind(M_trimmed_31_86_KOP_male_counts_only.txt, M_trimmed_31_86_KOP_queen_counts_only.txt)
m3189kop<-rbind(M_trimmed_31_89_KOP_male_counts_only.txt, M_trimmed_31_89_KOP_queen_counts_only.txt)


q0210kop<-rbind(Q_trimmed_02_10_KOP_male_counts_only.txt, Q_trimmed_02_10_KOP_queen_counts_only.txt)
q0212kop<-rbind(Q_trimmed_02_12_KOP_male_counts_only.txt, Q_trimmed_02_12_KOP_queen_counts_only.txt)
q0218kop<-rbind(Q_trimmed_02_18_KOP_male_counts_only.txt, Q_trimmed_02_18_KOP_queen_counts_only.txt)
q0222kop<-rbind(Q_trimmed_02_22_KOP_male_counts_only.txt, Q_trimmed_02_22_KOP_queen_counts_only.txt)
q0235kop<-rbind(Q_trimmed_02_35_KOP_male_counts_only.txt, Q_trimmed_02_35_KOP_queen_counts_only.txt)
q0237kop<-rbind(Q_trimmed_02_37_KOP_male_counts_only.txt, Q_trimmed_02_37_KOP_queen_counts_only.txt)
q0244kop<-rbind(Q_trimmed_02_44_KOP_male_counts_only.txt, Q_trimmed_02_44_KOP_queen_counts_only.txt)
q0254kop<-rbind(Q_trimmed_02_54_KOP_male_counts_only.txt, Q_trimmed_02_54_KOP_queen_counts_only.txt)
q1274kop<-rbind(Q_trimmed_12_74_KOP_male_counts_only.txt, Q_trimmed_12_74_KOP_queen_counts_only.txt)
q1275kop<-rbind(Q_trimmed_12_75_KOP_male_counts_only.txt, Q_trimmed_12_75_KOP_queen_counts_only.txt)
q1276kop<-rbind(Q_trimmed_12_76_KOP_male_counts_only.txt, Q_trimmed_12_76_KOP_queen_counts_only.txt)
q1278kop<-rbind(Q_trimmed_12_78_KOP_male_counts_only.txt, Q_trimmed_12_78_KOP_queen_counts_only.txt)
q1282kop<-rbind(Q_trimmed_12_82_KOP_male_counts_only.txt, Q_trimmed_12_82_KOP_queen_counts_only.txt)
q1284kop<-rbind(Q_trimmed_12_84_KOP_male_counts_only.txt, Q_trimmed_12_84_KOP_queen_counts_only.txt)
q1285kop<-rbind(Q_trimmed_12_85_KOP_male_counts_only.txt, Q_trimmed_12_85_KOP_queen_counts_only.txt)
q1289kop<-rbind(Q_trimmed_12_89_KOP_male_counts_only.txt, Q_trimmed_12_89_KOP_queen_counts_only.txt)
q2209kop<-rbind(Q_trimmed_22_09_KOP_male_counts_only.txt, Q_trimmed_22_09_KOP_queen_counts_only.txt)
q2218kop<-rbind(Q_trimmed_22_18_KOP_male_counts_only.txt, Q_trimmed_22_18_KOP_queen_counts_only.txt)
q2221kop<-rbind(Q_trimmed_22_21_KOP_male_counts_only.txt, Q_trimmed_22_21_KOP_queen_counts_only.txt)
q2223kop<-rbind(Q_trimmed_22_23_KOP_male_counts_only.txt, Q_trimmed_22_23_KOP_queen_counts_only.txt)
q2226kop<-rbind(Q_trimmed_22_26_KOP_male_counts_only.txt, Q_trimmed_22_26_KOP_queen_counts_only.txt)
q2230kop<-rbind(Q_trimmed_22_30_KOP_male_counts_only.txt, Q_trimmed_22_30_KOP_queen_counts_only.txt)
q2235kop<-rbind(Q_trimmed_22_35_KOP_male_counts_only.txt, Q_trimmed_22_35_KOP_queen_counts_only.txt)
q2249kop<-rbind(Q_trimmed_22_49_KOP_male_counts_only.txt, Q_trimmed_22_49_KOP_queen_counts_only.txt)
q3144kop<-rbind(Q_trimmed_31_44_KOP_male_counts_only.txt, Q_trimmed_31_44_KOP_queen_counts_only.txt)
q3145kop<-rbind(Q_trimmed_31_45_KOP_male_counts_only.txt, Q_trimmed_31_45_KOP_queen_counts_only.txt)
q3164kop<-rbind(Q_trimmed_31_64_KOP_male_counts_only.txt, Q_trimmed_31_64_KOP_queen_counts_only.txt)
q3165kop<-rbind(Q_trimmed_31_65_KOP_male_counts_only.txt, Q_trimmed_31_65_KOP_queen_counts_only.txt)
q3179kop<-rbind(Q_trimmed_31_79_KOP_male_counts_only.txt, Q_trimmed_31_79_KOP_queen_counts_only.txt)
q3184kop<-rbind(Q_trimmed_31_84_KOP_male_counts_only.txt, Q_trimmed_31_84_KOP_queen_counts_only.txt)
q3186kop<-rbind(Q_trimmed_31_86_KOP_male_counts_only.txt, Q_trimmed_31_86_KOP_queen_counts_only.txt)
q3189kop<-rbind(Q_trimmed_31_89_KOP_male_counts_only.txt, Q_trimmed_31_89_KOP_queen_counts_only.txt)

q0210abd<-rbind(Q_trimmed_02_10_ABD_male_counts_only.txt, Q_trimmed_02_10_ABD_queen_counts_only.txt)
q0212abd<-rbind(Q_trimmed_02_12_ABD_male_counts_only.txt, Q_trimmed_02_12_ABD_queen_counts_only.txt)
q0218abd<-rbind(Q_trimmed_02_18_ABD_male_counts_only.txt, Q_trimmed_02_18_ABD_queen_counts_only.txt)
q0222abd<-rbind(Q_trimmed_02_22_ABD_male_counts_only.txt, Q_trimmed_02_22_ABD_queen_counts_only.txt)
q0235abd<-rbind(Q_trimmed_02_35_ABD_male_counts_only.txt, Q_trimmed_02_35_ABD_queen_counts_only.txt)
q0237abd<-rbind(Q_trimmed_02_37_ABD_male_counts_only.txt, Q_trimmed_02_37_ABD_queen_counts_only.txt)
q0244abd<-rbind(Q_trimmed_02_44_ABD_male_counts_only.txt, Q_trimmed_02_44_ABD_queen_counts_only.txt)
q0254abd<-rbind(Q_trimmed_02_54_ABD_male_counts_only.txt, Q_trimmed_02_54_ABD_queen_counts_only.txt)
q1274abd<-rbind(Q_trimmed_12_74_ABD_male_counts_only.txt, Q_trimmed_12_74_ABD_queen_counts_only.txt)
q1275abd<-rbind(Q_trimmed_12_75_ABD_male_counts_only.txt, Q_trimmed_12_75_ABD_queen_counts_only.txt)
q1276abd<-rbind(Q_trimmed_12_76_ABD_male_counts_only.txt, Q_trimmed_12_76_ABD_queen_counts_only.txt)
q1278abd<-rbind(Q_trimmed_12_78_ABD_male_counts_only.txt, Q_trimmed_12_78_ABD_queen_counts_only.txt)
q1282abd<-rbind(Q_trimmed_12_82_ABD_male_counts_only.txt, Q_trimmed_12_82_ABD_queen_counts_only.txt)
q1284abd<-rbind(Q_trimmed_12_84_ABD_male_counts_only.txt, Q_trimmed_12_84_ABD_queen_counts_only.txt)
q1285abd<-rbind(Q_trimmed_12_85_ABD_male_counts_only.txt, Q_trimmed_12_85_ABD_queen_counts_only.txt)
q1289abd<-rbind(Q_trimmed_12_89_ABD_male_counts_only.txt, Q_trimmed_12_89_ABD_queen_counts_only.txt)
q2209abd<-rbind(Q_trimmed_22_09_ABD_male_counts_only.txt, Q_trimmed_22_09_ABD_queen_counts_only.txt)
q2218abd<-rbind(Q_trimmed_22_18_ABD_male_counts_only.txt, Q_trimmed_22_18_ABD_queen_counts_only.txt)
q2221abd<-rbind(Q_trimmed_22_21_ABD_male_counts_only.txt, Q_trimmed_22_21_ABD_queen_counts_only.txt)
q2223abd<-rbind(Q_trimmed_22_23_ABD_male_counts_only.txt, Q_trimmed_22_23_ABD_queen_counts_only.txt)
q2226abd<-rbind(Q_trimmed_22_26_ABD_male_counts_only.txt, Q_trimmed_22_26_ABD_queen_counts_only.txt)
q2230abd<-rbind(Q_trimmed_22_30_ABD_male_counts_only.txt, Q_trimmed_22_30_ABD_queen_counts_only.txt)
q2235abd<-rbind(Q_trimmed_22_35_ABD_male_counts_only.txt, Q_trimmed_22_35_ABD_queen_counts_only.txt)
q2249abd<-rbind(Q_trimmed_22_49_ABD_male_counts_only.txt, Q_trimmed_22_49_ABD_queen_counts_only.txt)
q3144abd<-rbind(Q_trimmed_31_44_ABD_male_counts_only.txt, Q_trimmed_31_44_ABD_queen_counts_only.txt)
q3145abd<-rbind(Q_trimmed_31_45_ABD_male_counts_only.txt, Q_trimmed_31_45_ABD_queen_counts_only.txt)
q3164abd<-rbind(Q_trimmed_31_64_ABD_male_counts_only.txt, Q_trimmed_31_64_ABD_queen_counts_only.txt)
q3165abd<-rbind(Q_trimmed_31_65_ABD_male_counts_only.txt, Q_trimmed_31_65_ABD_queen_counts_only.txt)
q3179abd<-rbind(Q_trimmed_31_79_ABD_male_counts_only.txt, Q_trimmed_31_79_ABD_queen_counts_only.txt)
q3184abd<-rbind(Q_trimmed_31_84_ABD_male_counts_only.txt, Q_trimmed_31_84_ABD_queen_counts_only.txt)
q3186abd<-rbind(Q_trimmed_31_86_ABD_male_counts_only.txt, Q_trimmed_31_86_ABD_queen_counts_only.txt)
q3189abd<-rbind(Q_trimmed_31_89_ABD_male_counts_only.txt, Q_trimmed_31_89_ABD_queen_counts_only.txt)

# Couldn't get a for loop to work as sql wouldn't take 'i' as a variable name, so changed the name
# manually for now, changing 'maternal_count' to paternal when needed etc.

x <- sqldf("SELECT sg.V2,
              sg.V3,
              fg.V4
              FROM q3189abd AS sg
              LEFT JOIN final_genes_sql AS fg 
              ON sg.V1 = fg.V1
              AND (sg.V2 >= fg.V2 AND sg.V2 <= fg.V3)")
x<-subset(x, !x$V4=="NA")
colnames(x)<-c("snp","maternal_count","geneID")
write.table(x, file="q3189abd_withGene.txt", sep="\t",row.names = F,quote=F)







