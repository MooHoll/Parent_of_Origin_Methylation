#!/bin/bash

#PBS -N alignment_to_ref
#PBS -l walltime=45:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load bowtie2/2.3.5.1
module load samtools/1.9

# Before alignment, make sure the reference genome has been processed by bismark:
# /scratch/monoallelic/hjm31/bin/Bismark-0.22.3/bismark_genome_preparation /scratch/monoallelic/hjm31/honeybee_ASM/genomes/honeybee

for file in $(ls *_1.fastq)
do
  	base=$(basename $file "_1.fastq")
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark --multicore 6 \
    /scratch/monoallelic/hjm32/asm_honeybee/genomes/honeybee \
     -1 ${base}_1.fastq -2 ${base}_2.fastq
done