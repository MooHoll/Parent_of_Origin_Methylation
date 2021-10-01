# Make a folder for each new reference genome and place .fa file inside

#!/bin/bash

#PBS -N genome_prep
#PBS -l walltime=02:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=1

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load bowtie2/2.3.5.1 
module load samtools/1.9

# put each N masked genome into it's own folder first
for folder in $(ls -d *_genome)
do
	/scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_genome_preparation ./${folder}
done

# -----------------------------------
# Alignment to N masked genomes
# -----------------------------------

#!/bin/bash

#PBS -N alignment
#PBS -l walltime=48:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=20

# Run script in the working directory it was submitted in
cd $PBS_O_WORKDIR 

# Load software needed
module load bowtie2/2.3.5.1 
module load samtools/1.9

for file in $(ls *1.fq.gz )
do
    base=$(basename ${file} "_trim_1.fq.gz")
    genome=${base:0:3}
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark --multicore 3 --non_bs_mm \
    /scratch/monoallelic/hjm32/honeybee_reanalysis/genomes/honeybee/N_masked/${genome}_genome \
    -1 ${base}_trim_1.fq.gz -2 ${base}_trim_2.fq.gz
done

 
for file in $(ls *.bam)
do
    /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/deduplicate_bismark -p ${file}
done
 
 
 
 
 
 