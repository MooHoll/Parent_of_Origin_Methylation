# Make a folder for each new reference genome and place .fa file inside

module load bowtie2/2.3.5.1 
module load samtools/1.9

for folder in $(ls -d *_genome)
do
	/scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_genome_preparation ./${folder}
done

# -----------------------------------
# Alignment to N masked genomes
# -----------------------------------

#!/bin/bash

#PBS -N alignment
#PBS -l walltime=08:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=20

# Run script in the working directory it was submitted in
cd $PBS_O_WORKDIR 

# Load software needed
module load bowtie2/2.3.5.1 
module load samtools/1.9

/scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark --multicore 3 --non_bs_mm \
/scratch/monoallelic/hjm32/bumblebee/genome/N_masked/08_genome \
 -1 trim_w08_1.fq.gz -2 trim_w08_2.fq.gz

/scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark --multicore 3 --non_bs_mm \
/scratch/monoallelic/hjm32/bumblebee/genome/N_masked/19_genome \
 -1 trim_w19_1.fq.gz -2 trim_w19_2.fq.gz

/scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark --multicore 3 --non_bs_mm \
/scratch/monoallelic/hjm32/bumblebee/genome/N_masked/23_genome \
 -1 trim_w23_1.fq.gz -2 trim_w23_2.fq.gz

/scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark --multicore 3 --non_bs_mm \
/scratch/monoallelic/hjm32/bumblebee/genome/N_masked/37_genome \
 -1 trim_w37_1.fq.gz -2 trim_w37_2.fq.gz

 
for file in $(ls *.bam)
do
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/deduplicate_bismark -p ${file}
done
 
 
 
 
 
 