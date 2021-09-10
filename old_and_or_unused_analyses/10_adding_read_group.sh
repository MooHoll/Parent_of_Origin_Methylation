#!/bin/bash

#PBS -N sorting_bams
#PBS -l walltime=02:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load samtools/1.3.2
module load java/1.8
module load picard/2.6.0

# Add read-group information to the alignments (needed for SNP calling)
for file in $(ls *.bam)
do
  	base=$(basename ${file} "_trimmed_deduplicated_sorted.bam")
    java -jar /cm/shared/apps/picard/2.6.0/picard.jar \
    AddOrReplaceReadGroups \
    I=${file} \
    O=${base}_trim_dedup_sort_RG.bam \
    RGID=0001${base} \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=NA \
    RGSM=${base}
done

## NOTE: can just run on the login nodes, ~7min per sample
## Also need to index and sort all files, can do on login node (may need an interactive job if memory too high)
#samtools index <file.bam>  
