# Subset the original alignments to the reference genome to only keep reads which have
# an informative SNP, i.e. a SNP unique to either the mother or the father

module load bedtools/2.28.0
module load samtools/1.9

# Count original number of aligned reads
samtools view -c -F 260 w08_trimmed_deduplicated_sorted.bam # 24258152
samtools view -c -F 260 w19_trimmed_deduplicated_sorted.bam # 25304366
samtools view -c -F 260 w23_trimmed_deduplicated_sorted.bam # 25825600
samtools view -c -F 260 w37_trimmed_deduplicated_sorted.bam # 24754058

# Take reads with SNPs only found in the fathers
bedtools intersect -wa -abam w08_trimmed_deduplicated_sorted.bam -b m08_final.vcf > w08_male_reads.bam
samtools view -c -F 260 w08_male_reads.bam # 3060715

bedtools intersect -wa -abam w19_trimmed_deduplicated_sorted.bam -b m19_final.vcf > w19_male_reads.bam
samtools view -c -F 260 w19_male_reads.bam # 3068164

bedtools intersect -wa -abam w23_trimmed_deduplicated_sorted.bam -b m23_final.vcf > w23_male_reads.bam
samtools view -c -F 260 w23_male_reads.bam # 3191733

bedtools intersect -wa -abam w37_trimmed_deduplicated_sorted.bam -b m37_final.vcf > w37_male_reads.bam
samtools view -c -F 260 w37_male_reads.bam # 3311999

# Take reads with SNPs only found in the mothers
bedtools intersect -wa -abam w08_trimmed_deduplicated_sorted.bam -b q08_final.vcf > w08_queen_reads.bam
samtools view -c -F 260 w08_queen_reads.bam # 

bedtools intersect -wa -abam w19_trimmed_deduplicated_sorted.bam -b q19_final.vcf > w19_queen_reads.bam
samtools view -c -F 260 w19_queen_reads.bam #

bedtools intersect -wa -abam w23_trimmed_deduplicated_sorted.bam -b q23_final.vcf > w23_queen_reads.bam
samtools view -c -F 260 w23_queen_reads.bam #

bedtools intersect -wa -abam w37_trimmed_deduplicated_sorted.bam -b q37_final.vcf > w37_queen_reads.bam
samtools view -c -F 260 w37_queen_reads.bam #