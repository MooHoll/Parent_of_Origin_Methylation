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

# Use single ended mode because the reads are only maternal/paternal of origin and 
# it's therefore no longer paired (usually you would use -p for paired-end mode)
for file in $(ls *_reads.bam)
do
    base=$(basename ${file} "_reads.bam")
    genome=${base:1:2}
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor \
    -s --no_overlap --comprehensive --bedgraph --report --cytosine_report \
    --genome_folder /scratch/monoallelic/hjm32/bumblebee/genome/old_genome/N_masked/${genome}_genome \
    ${file}
done

for file in $(ls *cov.gz)
do
    base=$(basename ${file} "_reads.bismark.cov.gz")
    genome=${base:1:2}
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/coverage2cytosine \
    -o ${base} --merge_CpGs \
    --genome_folder /scratch/monoallelic/hjm32/bumblebee/genome/old_genome/N_masked/${genome}_genome \
    ${file}
done

# To read into methylkit:
#CPGRaw <- methRead(sample.list, 
#                   sample.id = list("sample1", "sample2","sample3","sample4"),
#                   assembly="bter_1.0",
#                   treatment=c(0,0,1,1),
#                   context="CpG",
#                   dbtype = NA,
#                   pipeline = "bismarkCoverage",
#                   header = T, 
#                   sep = "\t", 
#                   mincov=1,
#                   dbdir= path)