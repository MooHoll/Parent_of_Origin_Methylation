#!/bin/bash

#PBS -N trimming
#PBS -l walltime=00:05:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16
#PBS -q devel

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR

# Load software needed (pigz needed for parallel zipping)
module load cutadapt/1.11

# Trim 10 bases from reads 
for file in $(ls *1.fastq)
do
	base=$(basename $file "_1.fastq")
	cutadapt \
	-u 10 -U 10 \
	-o trim_${base}_1.fastq\
	-p trim_${base}_2.fastq \
	${base}_1.fastq \
	${base}_2.fastq
done
