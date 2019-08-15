#!/bin/bash

#PBS -N m_extraction
#PBS -l walltime=30:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run in current working directory
cd $PBS_O_WORKDIR

# Load modules
module load bismark/0.18.1
module load samtools/1.3.2


# Extract context specific methylation information and remove tail end of read 2 after viewing 
# methylation bias plots. Keep the .bedgraph, .cov and C-context reports from this step 

for file in $(ls *.bam)
do
    bismark_methylation_extractor -p --comprehensive --bedgraph --genome_folder /scratch/monoallelic/hm257/genome/ ${file}
done


## Try again with --cutoff of 3?
## The minimum number of times any methylation state (methylated or unmethylated) has to be seen for a nucleotide before its methylation percentage is reported. Default: 1.
