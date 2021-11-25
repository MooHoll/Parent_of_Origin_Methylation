#!/bin/bash

#PBS -N fastqc
#PBS -l walltime=18:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Change directory to the one the job was submitted in
cd $PBS_O_WORKDIR 

# Load required modules
module load fastqc/0.11.5                                                                

# Run fastqc                                                                                                          
for file in $(ls *.fastq)
do
	fastqc -t 8 ${file}
done