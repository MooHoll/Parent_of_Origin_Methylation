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
module load bismark/0.18.1

for file in $(ls *.sam)
do
    deduplicate_bismark -p ${file}
done

grep 'XB:Z:0\|XA:Z:0\|@' male.trim_w08_1_bismark_bt2_pe.deduplicated.sam > m08_no_mismatched.sam
grep 'XB:Z:0\|XA:Z:0\|@' queen.trim_w08_1_bismark_bt2_pe.deduplicated.sam > q08_no_mismatched.sam

grep 'XB:Z:0\|XA:Z:0\|@' male.trim_w19_1_bismark_bt2_pe.deduplicated.sam > m19_no_mismatched.sam
grep 'XB:Z:0\|XA:Z:0\|@' queen.trim_w19_1_bismark_bt2_pe.deduplicated.sam > q19_no_mismatched.sam

grep 'XB:Z:0\|XA:Z:0\|@' male.trim_w23_1_bismark_bt2_pe.deduplicated.sam > m23_no_mismatched.sam
grep 'XB:Z:0\|XA:Z:0\|@' queen.trim_w23_1_bismark_bt2_pe.deduplicated.sam > q23_no_mismatched.sam

grep 'XB:Z:0\|XA:Z:0\|@' male.trim_w37_1_bismark_bt2_pe.deduplicated.sam > m37_no_mismatched.sam
grep 'XB:Z:0\|XA:Z:0\|@' queen.trim_w37_1_bismark_bt2_pe.deduplicated.sam > q37_no_mismatched.sam

# Outputs of the counting is done in the .o job file
for file in $(ls *pe.deduplicated.sam)
do
    samtools view -F 0x4 ${file} | cut -f1 | sort | uniq | wc -l
done

#male.trim_w08_1_bismark_bt2_pe.deduplicated.sam 12131230
#male.trim_w19_1_bismark_bt2_pe.deduplicated.sam 12653928
#male.trim_w23_1_bismark_bt2_pe.deduplicated.sam 12914084
#male.trim_w37_1_bismark_bt2_pe.deduplicated.sam 12378420
#queen.trim_w08_1_bismark_bt2_pe.deduplicated.sam 12129818
#queen.trim_w19_1_bismark_bt2_pe.deduplicated.sam 12654402
#queen.trim_w23_1_bismark_bt2_pe.deduplicated.sam 12914285
#queen.trim_w37_1_bismark_bt2_pe.deduplicated.sam 12378003

for file in $(ls *no_mismatches.sam)
do
    samtools view -F 0x4 ${file} | cut -f1 | sort | uniq | wc -l
done



###------------------------------------------------------------------
# Methylkit told me this is the right way to sort a sam 
# This didn't work, the headers got mixed into the reads
#for file in $(ls *no_mismatches.sam)
#do
#  	base=$(basename ${file} "no_mismatches.sam")
#    grep -v \'^[[:space:]]*\@\' ${file} \
#    | sort -k3,3 -k4,4n > ${base}no_mismatches_deduplicated_sorted.sam
#done

#m08_no_mismatches.sam 9436654
#m19_no_mismatches.sam 9949163
#m23_no_mismatches.sam 10215217
#m37_no_mismatches.sam 9764004
#q08_no_mismatches.sam 9386969
#q19_no_mismatches.sam 9948178
#q23_no_mismatches.sam 10202390
#q37_no_mismatches.sam 9731659

## DO BELOW FOR ALL FLILES
# Get the original header from the unsorted files
samtools view -H m08__no_mismatches_deduplicated_sorted.sam > header_m08.sam

# Remove all header info from the messed up files (checked by wc -l this only removes headers nothing else)
grep -v "@" m08_no_mismatches_deduplicated_sorted.sam > m08_noheaderinfo.sam

# Re-add header information to the sorted cleaned files (Can't use samtools reheader, only work to reheader bam)
module load picard/2.1.0 

for file in $(ls *_noheaderinfo.sam)
do
    base=$(basename ${file} "_noheaderinfo.sam")
    java -jar /cm/shared/apps/picard/2.1.0/picard.jar ReplaceSamHeader \
    I= ${file} HEADER= header_${base}.sam O= ${base}_no_mismatches_deduplicated_sorted.sam
done






