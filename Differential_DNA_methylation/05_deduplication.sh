#!/bin/bash

#PBS -N deduplicating
#PBS -l walltime=06:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in (80hrs)
cd $PBS_O_WORKDIR

# Load software needed
module load samtools/1.9
module load bowtie2/2.3.5.1

REF_FA=/scratch/monoallelic/hjm32/bumblebee/genome/old_genome

# Dedupliate all bam files, also gives a .txt report with how much data was removed
for file in $(ls *bam)
do
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/deduplicate_bismark -p ${file}
done

# make m-bias plots
for file in $(ls *deduplicated.bam)
do
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor -p \
    --mbias_only --report ${file}
done

# Generate bismark downstream files
for file in $(ls *deduplicated.bam)
do
  	base=$(basename ${file} ".deduplicated.bam")
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor \
    -p --no_overlap --comprehensive --bedgraph --report --cytosine_report \
    --genome_folder ${REF_FA} \
    ${file}
done

# Make destranded file for methylkit
for file in $(ls *cov.gz)
do
    base=$(basename ${file} ".bismark.cov.gz")
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/coverage2cytosine \
    -o ${base} --merge_CpGs \
    --genome_folder ${REF_FA} \
    ${file}
done