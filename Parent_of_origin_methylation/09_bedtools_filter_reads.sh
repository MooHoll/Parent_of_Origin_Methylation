# Subset the N masked new genome alignments to only keep reads which have
# an informative SNP, i.e. a SNP unique to either the mother or the father

module load bedtools/2.28.0
module load samtools/1.9

# Count original number of aligned reads
samtools view -c -F 260 trim_w08_1_bismark_bt2_pe.deduplicated.bam # 24318976
samtools view -c -F 260 trim_w19_1_bismark_bt2_pe.deduplicated.bam # 25238424
samtools view -c -F 260 trim_w23_1_bismark_bt2_pe.deduplicated.bam # 25865874
samtools view -c -F 260 trim_w37_1_bismark_bt2_pe.deduplicated.bam # 24794764

# Need to unzip all vcfs for this to work, urgh
module load htslib/1.9
bgzip -d m08_final.vcf.gz
bgzip -d m19_final.vcf.gz
bgzip -d m23_final.vcf.gz
bgzip -d m37_final.vcf.gz
bgzip -d q08_final.vcf.gz
bgzip -d q19_final.vcf.gz
bgzip -d q23_final.vcf.gz
bgzip -d q37_final.vcf.gz

# Take reads with SNPs only found in the fathers
bedtools intersect -wa -abam trim_w08_1_bismark_bt2_pe.deduplicated.bam -b ../snps/m08_final.vcf > w08_male_reads.bam
samtools view -c -F 260 w08_male_reads.bam # 3207495 (3057385) new genome

bedtools intersect -wa -abam trim_w19_1_bismark_bt2_pe.deduplicated.bam -b ../snps/m19_final.vcf > w19_male_reads.bam
samtools view -c -F 260 w19_male_reads.bam # 3232624 (3138721)

bedtools intersect -wa -abam trim_w23_1_bismark_bt2_pe.deduplicated.bam -b ../snps/m23_final.vcf > w23_male_reads.bam
samtools view -c -F 260 w23_male_reads.bam # 3337570 (3100481)

bedtools intersect -wa -abam trim_w37_1_bismark_bt2_pe.deduplicated.bam -b ../snps/m37_final.vcf > w37_male_reads.bam
samtools view -c -F 260 w37_male_reads.bam # 3459605 (3393155)

# Take reads with SNPs only found in the mothers
bedtools intersect -wa -abam trim_w08_1_bismark_bt2_pe.deduplicated.bam -b ../snps/q08_final.vcf > w08_queen_reads.bam
samtools view -c -F 260 w08_queen_reads.bam # 960798 (1049566)

bedtools intersect -wa -abam trim_w19_1_bismark_bt2_pe.deduplicated.bam -b ../snps/q19_final.vcf > w19_queen_reads.bam
samtools view -c -F 260 w19_queen_reads.bam # 922541 (948770)

bedtools intersect -wa -abam trim_w23_1_bismark_bt2_pe.deduplicated.bam -b ../snps/q23_final.vcf > w23_queen_reads.bam
samtools view -c -F 260 w23_queen_reads.bam # 949697 (1073333)

bedtools intersect -wa -abam trim_w37_1_bismark_bt2_pe.deduplicated.bam -b ../snps/q37_final.vcf > w37_queen_reads.bam
samtools view -c -F 260 w37_queen_reads.bam # 872826 (873228)

for file in $(ls *reads.bam)
do
  	base=$(basename $file "_reads.bam")
    samtools sort -o ${base}_reads_sorted.bam ${file}
done