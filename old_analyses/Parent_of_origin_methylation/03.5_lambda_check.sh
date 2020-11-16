#!/bin/bash

#PBS -N alignment_to_lambda
#PBS -l walltime=04:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in
cd $PBS_O_WORKDIR

# Load software needed
module load bismark/0.18.1
module load bowtie2/2.2.9
module load samtools/1.3.2

for file in $(ls *1.fq.gz)
do
  	base=$(basename $file "1.fq.gz")
    bismark --multicore 3 /scratch/monoallelic/hm257/repro_methylation/merged/lambda_genome -1 ${base}1.fq.gz -2 ${base}2.fq.gz
done