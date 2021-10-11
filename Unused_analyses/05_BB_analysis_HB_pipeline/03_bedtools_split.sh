# Subset the alignments to the N masked genome to only keep reads which have
# an informative SNP, i.e. a SNP unique to either the mother or the father

module load bedtools/2.28.0
module load samtools/1.9

# Count original number of aligned reads
samtools view -c -F 260 trim_w08_1_bismark_bt2_pe.deduplicated.bam # 24666038
samtools view -c -F 260 trim_w19_1_bismark_bt2_pe.deduplicated.bam # 25718716
samtools view -c -F 260 trim_w23_1_bismark_bt2_pe.deduplicated.bam # 26223572
samtools view -c -F 260 trim_w37_1_bismark_bt2_pe.deduplicated.bam # 25151300

#NOTE: alignments slightl improved compared to normal and alternative ref genomes with SNP placed

# Take reads with SNPs only found in the fathers
bedtools intersect -wa -abam trim_w08_1_bismark_bt2_pe.deduplicated.bam -b ../snps/m08_final.vcf > w08_male_reads.bam
samtools view -c -F 260 w08_male_reads.bam # 3301564

bedtools intersect -wa -abam trim_w19_1_bismark_bt2_pe.deduplicated.bam -b ../snps/m19_final.vcf > w19_male_reads.bam
samtools view -c -F 260 w19_male_reads.bam # 3325469

bedtools intersect -wa -abam trim_w23_1_bismark_bt2_pe.deduplicated.bam -b ../snps/m23_final.vcf > w23_male_reads.bam
samtools view -c -F 260 w23_male_reads.bam # 3422310

bedtools intersect -wa -abam trim_w37_1_bismark_bt2_pe.deduplicated.bam -b ../snps/m37_final.vcf > w37_male_reads.bam
samtools view -c -F 260 w37_male_reads.bam # 3558016

# Take reads with SNPs only found in the mothers
bedtools intersect -wa -abam trim_w08_1_bismark_bt2_pe.deduplicated.bam -b ../snps/q08_final.vcf > w08_queen_reads.bam
samtools view -c -F 260 w08_queen_reads.bam # 982909

bedtools intersect -wa -abam trim_w19_1_bismark_bt2_pe.deduplicated.bam -b ../snps/q19_final.vcf > w19_queen_reads.bam
samtools view -c -F 260 w19_queen_reads.bam # 939155

bedtools intersect -wa -abam trim_w23_1_bismark_bt2_pe.deduplicated.bam -b ../snps/q23_final.vcf > w23_queen_reads.bam
samtools view -c -F 260 w23_queen_reads.bam # 969378

bedtools intersect -wa -abam trim_w37_1_bismark_bt2_pe.deduplicated.bam -b ../snps/q37_final.vcf > w37_queen_reads.bam
samtools view -c -F 260 w37_queen_reads.bam # 894990

for file in $(ls *reads.bam)
do
  	base=$(basename $file "_reads.bam")
    samtools sort -o ${base}_reads_sorted.bam ${file}
done