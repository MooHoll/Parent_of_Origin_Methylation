#!/bin/bash

#PBS -N mbias_only
#PBS -l walltime=09:00:00 
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in
cd $PBS_O_WORKDIR

# Load software needed
module load samtools/1.9


# Produce m-bias plots to determine filtering during methylation_extraction commmand
REF_FA=/scratch/monoallelic/hm257/meth_leuven/genome

for file in $(ls *.bam)
do
  	base=$(basename ${file} ".deduplicated.bam")
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor -p \
    --mbias_only --report ${base}.deduplicated.bam
done