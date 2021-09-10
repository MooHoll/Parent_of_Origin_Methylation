#!/bin/bash

#PBS -N deduplicating_bams
#PBS -l walltime=02:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in
cd $PBS_O_WORKDIR

# Load software needed
module load samtools/1.9

for file in $(ls *deduplicated.bam)
do
    base=$(basename ${file} "_1_bismark_bt2_pe.deduplicated.bam")
    samtools view -h -o ${base}.sam ${file}
done

grep 'XB:Z:0\|XA:Z:0\|@' male.trim_w08.sam > male_w08_no_mismatches.sam
grep 'XB:Z:0\|XA:Z:0\|@' queen.trim_w08.sam > queen_w08_no_mismatches.sam

grep 'XB:Z:0\|XA:Z:0\|@' male.trim_w19.sam > male_w19_no_mismatches.sam
grep 'XB:Z:0\|XA:Z:0\|@' queen.trim_w19.sam > queen_w19_no_mismatches.sam

grep 'XB:Z:0\|XA:Z:0\|@' male.trim_w23.sam > male_w23_no_mismatches.sam
grep 'XB:Z:0\|XA:Z:0\|@' queen.trim_w23.sam > queen_w23_no_mismatches.sam

grep 'XB:Z:0\|XA:Z:0\|@' male.trim_w37.sam > male_w37_no_mismatches.sam
grep 'XB:Z:0\|XA:Z:0\|@' queen.trim_w37.sam > queen_w37_no_mismatches.sam

# Outputs of the counting is done in the .o job file
for file in $(ls *.sam)
do
    samtools view -F 0x4 ${file} | cut -f1 | sort | uniq | wc -l
done

#male.trim_w08.sam 12206962
#male.trim_w19.sam 12737961
#male.trim_w23.sam 12992925
#male.trim_w37.sam 12455055
#male_w08_no_mismatches.sam 9492258
#male_w19_no_mismatches.sam 10017542
#male_w23_no_mismatches.sam 10307951
#male_w37_no_mismatches.sam 9849215

#queen.trim_w08.sam 12166873
#queen.trim_w19.sam 12687077
#queen.trim_w23.sam 12950740
#queen.trim_w37.sam 12409875
#queen_w08_no_mismatches.sam 9365839
#queen_w19_no_mismatches.sam 9896961
#queen_w23_no_mismatches.sam 10174728
#queen_w37_no_mismatches.sam 9706917

#----------------------------------------------------------------------------------
## DO BELOW FOR ALL FLILES

# Methylkit told me this is the right way to sort a sam 

for file in $(ls *no_mismatches.sam)
do
  	base=$(basename ${file} "no_mismatches.sam")
    grep -v \'^[[:space:]]*\@\' ${file} \
    | sort -k3,3 -k4,4n > ${base}no_mismatches_sorted.sam
done

# It fucks it up though so need to fix it:

# Get the original header from the unsorted files
for file in $(ls *mismatches.sam)
do
    base=$(basename ${file} "_no_mismatches.sam")
    samtools view -H male_w08_no_mismatches.sam > header_${base}.sam
done

# Remove all header info from the messed up files (checked by wc -l this only removes headers nothing else)
for file in $(ls *no_mismatches_sorted.sam)
do
    base=$(basename ${file} "_no_mismatches_sorted.sam")
    grep -v "@" male_w08_no_mismatches_sorted.sam > ${base}_noheaderinfo.sam
done

# Re-add header information to the sorted cleaned files (Can't use samtools reheader, only work to reheader bam)
module load picard/2.1.0 

for file in $(ls *_noheaderinfo.sam)
do
    base=$(basename ${file} "_noheaderinfo.sam")
    java -jar /cm/shared/apps/picard/2.1.0/picard.jar ReplaceSamHeader \
    I=${file} HEADER=header_${base}.sam O=${base}_for_methylkit.sam
done






