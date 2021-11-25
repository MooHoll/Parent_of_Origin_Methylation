#!/bin/bash

#PBS -N getting_data
#PBS -l walltime=16:30:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Change directory to the one the job was submitted in
cd $PBS_O_WORKDIR 

# Load required modules
module load sratoolkit/2.11.1

prefetch --option-file accessions.txt

cat accessions.txt | while read line
do
fasterq-dump --split-files ${line}.sra
done