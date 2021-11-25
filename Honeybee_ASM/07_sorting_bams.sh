#!/bin/bash

#PBS -N sorting_bams
#PBS -l walltime=12:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load samtools/1.9

for file in $(ls *deduplicated.bam)
do
    base=$(basename $file "_1_bismark_bt2_pe.deduplicated.bam")
    samtools sort -n -@ 15 -O bam -T _tmp \
    -o ${base}_sorted_methtuple.bam ${file}
done