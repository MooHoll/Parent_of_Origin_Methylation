#!/bin/bash

#PBS -N alignment_to_lambda
#PBS -l walltime=06:00:00
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
# /scratch/monoallelic/hjm31/bin/Bismark-0.22.3/bismark_genome_preparation /scratch/monoallelic/hjm31/honeybee_ASM/genomes/lambda_genome

for file in $(ls *_1.fastq)
do
  	base=$(basename ${file} "_1.fastq")
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark --multicore 6 \
    --prefix lambda \
    /scratch/monoallelic/hjm32/asm_honeybee/genomes/lambda_genome \
     -1 ${base}_1.fastq -2 ${base}_2.fastq
done

cat to_do.txt | while read line
do
  	base=$(basename ${line} "_1.fastq")
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark --multicore 6 \
    --prefix lambda \
    /scratch/monoallelic/hjm32/asm_honeybee/genomes/lambda_genome \
     -1 ${base}_1.fastq -2 ${base}_2.fastq
done