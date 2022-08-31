# ------------------------------------------
# Count the number of SNPs per gene
# ------------------------------------------

# Make a gene bed file (chr, start, end, gene_id)
grep "gene" ref_Bter_1.0_top_level_numbered_exons.txt > genes.txt
awk 'BEGIN{OFS="\t"}{print $1,$3,$4,$6}' genes.txt > genes1.txt
sed '1d' genes1.txt > genes.bed

# Count all SNPs per gene for each sample
module load bedtools/2.28.0

for file in $(ls *filtered.recode.vcf)
do
    base=$(basename ${file} "_filtered.recode.vcf")
    bedtools intersect -a genes.bed \
    -b ${file} -c > ${base}_SNP_counts_per_gene.txt
done

# count number of snps total in output
awk '{s+=$5} END {print s}' m08_SNP_counts_per_gene.txt
# 493520
awk '{s+=$5} END {print s}' m19_SNP_counts_per_gene.txt
# 453586
awk '{s+=$5} END {print s}' m23_SNP_counts_per_gene.txt
# 446435
awk '{s+=$5} END {print s}' m37_SNP_counts_per_gene.txt
# 499655
awk '{s+=$5} END {print s}' q08_SNP_counts_per_gene.txt
# 695593
awk '{s+=$5} END {print s}' q19_SNP_counts_per_gene.txt
# 637714
awk '{s+=$5} END {print s}' q23_SNP_counts_per_gene.txt
# 632870
awk '{s+=$5} END {print s}' q37_SNP_counts_per_gene.txt
# 628635