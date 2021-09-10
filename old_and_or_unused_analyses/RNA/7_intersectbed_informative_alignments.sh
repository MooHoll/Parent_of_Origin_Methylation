#!/bin/bash

#PBS -N intersect_bed
#PBS -l walltime=03:30:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=2:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load bedtools/2.28.0

##############################
# Keep reads which have an informative SNP (i.e. a SNP unique to either the mother or father)


# To run a test:
#  qsub -A intersecting -I -l walltime=01:00:00 -l pvmem=10gb -l nodes=1:ppn=1

bedtools intersect -wa -abam m08_trimmed_deduplicated_sorted.bam -b m08_final.vcf -sorted > m08_reads_with_SNP.bam
# Total reads before
# Total reads after

