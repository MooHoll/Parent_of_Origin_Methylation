#!/bin/bash

#PBS -N alignment
#PBS -l walltime=00:05:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=20
#PBS -q devel

# Run script in the working directory it was submitted in (25hrs)
cd $PBS_O_WORKDIR 

# Load software needed
module load bowtie2/2.2.9
module load samtools/1.3.2

# Align all samples to the reference Bter_1.0
/scratch/monoallelic/hjm31/bin_new/Bismark-0.22.3/bismark_genome_preparation \
/scratch/monoallelic/hjm32/bumblebee/genome/new_genome

REF_FA=/scratch/monoallelic/hjm32/bumblebee/genome/new_genome

for file in $(ls *1.fq.gz)
do
	base=$(basename $file "1.fq.gz")
	/scratch/monoallelic/hjm31/bin_new/Bismark-0.22.3/bismark \
	--multicore 3 \
	${REF_FA} -1 ${base}1.fq.gz -2 ${base}2.fq.gz
done 