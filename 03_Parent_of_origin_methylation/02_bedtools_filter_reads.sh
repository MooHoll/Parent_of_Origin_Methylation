# Subset the alt_ref alignments to only keep reads which have
# an informative SNP, i.e. a SNP unique to either the mother or the father

module load bedtools/2.28.0
module load samtools/1.9

# Count original number of aligned reads
samtools view -c -F 260 male.trim_w08_1_bismark_bt2_pe.deduplicated.bam # 24413924
samtools view -c -F 260 male.trim_w19_1_bismark_bt2_pe.deduplicated.bam # 25475922
samtools view -c -F 260 male.trim_w23_1_bismark_bt2_pe.deduplicated.bam # 25985850
samtools view -c -F 260 male.trim_w37_1_bismark_bt2_pe.deduplicated.bam # 24910110
samtools view -c -F 260 queen.trim_w08_1_bismark_bt2_pe.deduplicated.bam # 24333746
samtools view -c -F 260 queen.trim_w19_1_bismark_bt2_pe.deduplicated.bam # 25374154
samtools view -c -F 260 queen.trim_w23_1_bismark_bt2_pe.deduplicated.bam # 25901480
samtools view -c -F 260 queen.trim_w37_1_bismark_bt2_pe.deduplicated.bam # 24819750

# Take reads with SNPs only found in the fathers
bedtools intersect -wa -abam male.trim_w08_1_bismark_bt2_pe.deduplicated.bam -b ../snps/m08_final.vcf > w08_male_reads.bam
samtools view -c -F 260 w08_male_reads.bam # 3207495

bedtools intersect -wa -abam male.trim_w19_1_bismark_bt2_pe.deduplicated.bam -b ../snps/m19_final.vcf > w19_male_reads.bam
samtools view -c -F 260 w19_male_reads.bam # 3232624

bedtools intersect -wa -abam male.trim_w23_1_bismark_bt2_pe.deduplicated.bam -b ../snps/m23_final.vcf > w23_male_reads.bam
samtools view -c -F 260 w23_male_reads.bam # 3337570

bedtools intersect -wa -abam male.trim_w37_1_bismark_bt2_pe.deduplicated.bam -b ../snps/m37_final.vcf > w37_male_reads.bam
samtools view -c -F 260 w37_male_reads.bam # 3459605

# Take reads with SNPs only found in the mothers
bedtools intersect -wa -abam queen.trim_w08_1_bismark_bt2_pe.deduplicated.bam -b ../snps/q08_final.vcf > w08_queen_reads.bam
samtools view -c -F 260 w08_queen_reads.bam # 960798

bedtools intersect -wa -abam queen.trim_w19_1_bismark_bt2_pe.deduplicated.bam -b ../snps/q19_final.vcf > w19_queen_reads.bam
samtools view -c -F 260 w19_queen_reads.bam # 922541

bedtools intersect -wa -abam queen.trim_w23_1_bismark_bt2_pe.deduplicated.bam -b ../snps/q23_final.vcf > w23_queen_reads.bam
samtools view -c -F 260 w23_queen_reads.bam # 949697

bedtools intersect -wa -abam queen.trim_w37_1_bismark_bt2_pe.deduplicated.bam -b ../snps/q37_final.vcf > w37_queen_reads.bam
samtools view -c -F 260 w37_queen_reads.bam # 872826