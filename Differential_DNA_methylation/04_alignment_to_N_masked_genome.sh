# Make a folder for each new reference genome and place .fa file inside

module load bowtie2/2.3.5.1 
module load samtools/1.9

for folder in $(ls -d *)
do
	/scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_genome_preparation ./${folder}
done

# -----------------------------------
# Alignment to N masked genomes
# -----------------------------------

#!/bin/bash

#PBS -N alignment
#PBS -l walltime=36:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=20

# Run script in the working directory it was submitted in
cd $PBS_O_WORKDIR 

# Load software needed
module load bowtie2/2.3.5.1 
module load samtools/1.9

for file in $(ls *_1.fq.gz)
do
    base=$(basename ${file} "_1.fq.gz")
    genome=${base:1:2}
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark --multicore 3 --non_bs_mm \
    /scratch/monoallelic/hjm32/bumblebee/masked_genomes/${genome}_genome \
    -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz
done