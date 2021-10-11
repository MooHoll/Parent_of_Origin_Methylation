# Subset the alignments to the N masked genome to only keep reads which have
# an informative SNP, i.e. a SNP unique to either the mother or the father

module load bedtools/2.28.0
module load samtools/1.9

# Count original number of aligned reads
for file in $(ls *.bam)
do
  samtools view -c -F 260 ${file}
done

#875_1_trim_1_bismark_bt2_pe.deduplicated.bam 64119602
#875_2_trim_1_bismark_bt2_pe.deduplicated.bam 67012318
#875_E_trim_1_bismark_bt2_pe.deduplicated.bam 62476834
#875_F_trim_1_bismark_bt2_pe.deduplicated.bam 67357370
#882_0_trim_1_bismark_bt2_pe.deduplicated.bam 48685362
#882_2_trim_1_bismark_bt2_pe.deduplicated.bam 83302064
#882_3_trim_1_bismark_bt2_pe.deduplicated.bam 86465966
#882_5_trim_1_bismark_bt2_pe.deduplicated.bam 56243164
#888_1_trim_1_bismark_bt2_pe.deduplicated.bam 33602012
#888_2_trim_1_bismark_bt2_pe.deduplicated.bam 43645122
#888_E_trim_1_bismark_bt2_pe.deduplicated.bam 70628674
#888_F_trim_1_bismark_bt2_pe.deduplicated.bam 48381844
#894_1_trim_1_bismark_bt2_pe.deduplicated.bam 77372098
#894_2_trim_1_bismark_bt2_pe.deduplicated.bam 83784522
#894_5_trim_1_bismark_bt2_pe.deduplicated.bam 84372392
#894_6_trim_1_bismark_bt2_pe.deduplicated.bam 76957098

# Upzip snp files
module load htslib/1.9
bgzip -d *final.vcf.gz

#-----------------------------------------------------------

#!/bin/bash

#PBS -N split_reads
#PBS -l walltime=10:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

module load bedtools/2.28.0
module load samtools/1.9

# Take reads with SNPs only found in the fathers
for file in $(ls *deduplicated.bam)
do
  base=$(basename ${file} "_trim_1_bismark_bt2_pe.deduplicated.bam")
  genome=${base:0:3}
  bedtools intersect -wa -abam ${file} -b ../../snps/${genome}_drone_final.vcf > ${base}_drone_reads.bam
  samtools view -c -F 260 ${base}_drone_reads.bam
done

#875_1_drone_reads.bam 12195828
#875_2_drone_reads.bam 13352993
#875_E_drone_reads.bam 12124381
#875_F_drone_reads.bam 13481486
#882_0_drone_reads.bam 8112766
#882_2_drone_reads.bam 12679041
#882_3_drone_reads.bam 14715829
#882_5_drone_reads.bam 8635016
#888_1_drone_reads.bam 2893803
#888_2_drone_reads.bam 4046505
#888_E_drone_reads.bam 6285500
#888_F_drone_reads.bam 3461273
#894_1_drone_reads.bam 8224916
#894_2_drone_reads.bam 7486348
#894_5_drone_reads.bam 10360638
#894_6_drone_reads.bam 9549647


# Take reads with SNPs only found in the mothers
for file in $(ls *deduplicated.bam)
do
  base=$(basename ${file} "_trim_1_bismark_bt2_pe.deduplicated.bam")
  genome=${base:0:3}
  bedtools intersect -wa -abam ${file} -b ../../snps/${genome}_queen_final.vcf > ${base}_queen_reads.bam
  samtools view -c -F 260 ${base}_queen_reads.bam
done

#875_1_queen_reads.bam 2360185
#875_2_queen_reads.bam 2517433
#875_E_queen_reads.bam 2259410
#875_F_queen_reads.bam 2450915
#882_0_queen_reads.bam 1521841
#882_2_queen_reads.bam 2314424
#882_3_queen_reads.bam 2801198
#882_5_queen_reads.bam 1574245
#888_1_queen_reads.bam 3296276
#888_2_queen_reads.bam 4438115
#888_E_queen_reads.bam 6703398
#888_F_queen_reads.bam 4442349
#894_1_queen_reads.bam 3687597
#894_2_queen_reads.bam 3062478
#894_5_queen_reads.bam 4587669
#894_6_queen_reads.bam 4292043


# Sort the files for methylkit
for file in $(ls *reads.bam)
do
  	base=$(basename $file "_reads.bam")
    samtools sort -@ 7 -o ${base}_reads_sorted.bam ${file}
done