#!/bin/bash

#PBS -N sorting_bams
#PBS -l walltime=02:30:00
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


# Sort all bams ready for SNP calling

for file in $(ls *deduplicated.bam)
do
  	base=$(basename $file "1_bismark_bt2_pe_deduplicated.bam")
    samtools sort -o ${base}deduplicated_sorted.bam ${file}
done
