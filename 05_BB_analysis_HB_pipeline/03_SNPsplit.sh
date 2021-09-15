#!/bin/bash

#PBS -N SNP_split
#PBS -l walltime=12:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=20

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load samtools/1.9

/scratch/monoallelic/hjm31/bin_new/SNPsplit-0.5.0/SNPsplit \
--paired --bisulfite \
--snp_file <X> \
<x>
 
for file in $(ls *.bam)
do
    /scratch/monoallelic/hjm31/bin_new/Bismark-0.22.3/deduplicate_bismark -p ${file}
done
 
 
 
 
 
 