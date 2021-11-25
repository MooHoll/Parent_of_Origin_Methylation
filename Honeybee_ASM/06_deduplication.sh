#!/bin/bash

#PBS -N deduplication
#PBS -l walltime=24:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load bowtie2/2.3.5.1
module load samtools/1.9

for file in $(ls *bam)
do
	/scratch/monoallelic/hjm32/bin/Bismark-0.22.3/deduplicate_bismark \
	-p --bam ${file}
done