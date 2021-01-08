# Make a folder for each new reference genome and place .fa file inside (takes a while)

module load bowtie2/2.3.5.1 
module load samtools/1.9
# Bismark-0.22.3

for folder in $(ls -d *_genome)
do
	/scratch/monoallelic/hm257/bb_imprinting/bin/Bismark-0.22.3/bismark_genome_preparation ./${folder}
done

# -----------------------------------
# Alignment to parental genomes
# -----------------------------------

#!/bin/bash

#PBS -N alignment_to_ref
#PBS -l walltime=23:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=20

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load bowtie2/2.3.5.1 
module load samtools/1.9


/scratch/monoallelic/hm257/bb_imprinting/bin/Bismark-0.22.3/bismark --multicore 3 --prefix male \
/scratch/monoallelic/hm257/bb_imprinting/alternate_references/m08_genome \
 -1 trim_w08_1.fq.gz -2 trim_w08_2.fq.gz
 
/scratch/monoallelic/hm257/bb_imprinting/bin/Bismark-0.22.3/bismark --multicore 3 --prefix queen \
/scratch/monoallelic/hm257/bb_imprinting/alternate_references/q08_genome \
 -1 trim_w08_1.fq.gz -2 trim_w08_2.fq.gz
 

 /scratch/monoallelic/hm257/bb_imprinting/bin/Bismark-0.22.3/bismark --multicore 3 --prefix male \
/scratch/monoallelic/hm257/bb_imprinting/alternate_references/m19_genome \
 -1 trim_w19_1.fq.gz -2 trim_w19_2.fq.gz
 
/scratch/monoallelic/hm257/bb_imprinting/bin/Bismark-0.22.3/bismark --multicore 3 --prefix queen \
/scratch/monoallelic/hm257/bb_imprinting/alternate_references/q19_genome \
 -1 trim_w19_1.fq.gz -2 trim_w19_2.fq.gz
 

/scratch/monoallelic/hm257/bb_imprinting/bin/Bismark-0.22.3/bismark --multicore 3 --prefix male \
/scratch/monoallelic/hm257/bb_imprinting/alternate_references/m23_genome \
 -1 trim_w23_1.fq.gz -2 trim_w23_2.fq.gz
 
/scratch/monoallelic/hm257/bb_imprinting/bin/Bismark-0.22.3/bismark --multicore 3 --prefix queen \
/scratch/monoallelic/hm257/bb_imprinting/alternate_references/q23_genome \
 -1 trim_w23_1.fq.gz -2 trim_w23_2.fq.gz

 
/scratch/monoallelic/hm257/bb_imprinting/bin/Bismark-0.22.3/bismark --multicore 3 --prefix male \
/scratch/monoallelic/hm257/bb_imprinting/alternate_references/m37_genome \
 -1 trim_w37_1.fq.gz -2 trim_w37_2.fq.gz
 
/scratch/monoallelic/hm257/bb_imprinting/bin/Bismark-0.22.3/bismark --multicore 3 --prefix queen \
/scratch/monoallelic/hm257/bb_imprinting/alternate_references/q37_genome \
 -1 trim_w37_1.fq.gz -2 trim_w37_2.fq.gz
 
 
 
 
 
 
 
 
 