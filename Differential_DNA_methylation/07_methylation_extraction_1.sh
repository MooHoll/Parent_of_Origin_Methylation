#!/bin/bash

#PBS -N m_extraction
#PBS -l walltime=16:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run in current working directory (30hrs)
cd $PBS_O_WORKDIR

# Load modules
module load samtools/1.9

REF_FA=/scratch/monoallelic/hjm32/bumblebee/genome/old_genome

for file in $(ls *.bam)
do
  	base=$(basename ${file} ".deduplicated.bam")
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor \
    -p --no_overlap --comprehensive --bedgraph --report --cytosine_report \
    --genome_folder ${REF_FA} \
    ${file}
done

for file in $(ls *cov.gz)
do
    base=$(basename ${file} ".bismark.cov.gz")
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/coverage2cytosine \
    -o ${base} --merge_CpGs \
    --genome_folder ${REF_FA} \
    ${file}
done