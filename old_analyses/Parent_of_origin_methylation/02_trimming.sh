#!/bin/bash

#PBS -N trimming
#PBS -l walltime=10:30:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed (pigz needed for parallel zipping)
module load cutadapt/1.11
module load pigz/2.3.3

# create the directory where the output files are to be written 
OUTPUT=trimmed_data                                                                                                                                     
if [ ! -d "$OUTPUT" ]; then                                                                                 
    mkdir -p ${OUTPUT} 
fi

# Trim 10 bases from both ends of both reads, also minimum read length must be 20bp
for file in $(ls *1.fq.gz)
do
	base=$(basename $file "1.fq.gz")
	cutadapt -u 10 -U 10 -m 20 -o ${OUTPUT}/trim_${base}1.fq.gz -p ${OUTPUT}/trim_${base}2.fq.gz ${base}1.fq.gz ${base}2.fq.gz
done 