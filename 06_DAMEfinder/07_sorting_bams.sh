#!/bin/bash

#PBS -N sorting_bams
#PBS -l walltime=02:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load samtools/1.9

# Using bams aligned to the Bter 1.0 ref genome
for file in $(ls *deduplicated_sorted.bam)
do
    base=$(basename $file "_trimmed_deduplicated_sorted.bam")
    samtools sort -n -@ 15 -O bam -T _tmp \
    -o ${base}_sorted_methtuple.bam ${file}
done