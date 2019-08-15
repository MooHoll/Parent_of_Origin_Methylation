# Make a folder for each new reference genome and place .fa file inside (takes a while)

module load bismark/0.18.1
module load bowtie2/2.2.9
module load samtools/1.3.2

for folder in $(ls -d *_genome)
do
	bismark_genome_preparation ./${folder}
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
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load bismark/0.18.1
module load bowtie2/2.2.9
module load samtools/1.3.2


# Write out additional information about the number of non-bisulfite mismatches
# to do this need SAM output format, which means can't have multicore option
bismark --non_bs_mm --sam --prefix male \
/scratch/monoallelic/hm257/leuven_meth/alternate_references/m08_genome \
 -1 trim_w08_1.fq.gz -2 trim_w08_2.fq.gz
 
bismark --non_bs_mm --sam --prefix queen \
/scratch/monoallelic/hm257/leuven_meth/alternate_references/q08_genome \
 -1 trim_w08_1.fq.gz -2 trim_w08_2.fq.gz
 
 
 
 bismark --non_bs_mm --sam --prefix male \
/scratch/monoallelic/hm257/leuven_meth/alternate_references/m19_genome \
 -1 trim_w19_1.fq.gz -2 trim_w19_2.fq.gz
 
bismark --non_bs_mm --sam --prefix queen \
/scratch/monoallelic/hm257/leuven_meth/alternate_references/q19_genome \
 -1 trim_w19_1.fq.gz -2 trim_w19_2.fq.gz
 
 
 
  bismark --non_bs_mm --sam --prefix male \
/scratch/monoallelic/hm257/leuven_meth/alternate_references/m23_genome \
 -1 trim_w23_1.fq.gz -2 trim_w23_2.fq.gz
 
bismark --non_bs_mm --sam --prefix queen \
/scratch/monoallelic/hm257/leuven_meth/alternate_references/q23_genome \
 -1 trim_w23_1.fq.gz -2 trim_w23_2.fq.gz
  
 
 
  bismark --non_bs_mm --sam --prefix male \
/scratch/monoallelic/hm257/leuven_meth/alternate_references/m37_genome \
 -1 trim_w37_1.fq.gz -2 trim_w37_2.fq.gz
 
bismark --non_bs_mm --sam --prefix queen \
/scratch/monoallelic/hm257/leuven_meth/alternate_references/q37_genome \
 -1 trim_w37_1.fq.gz -2 trim_w37_2.fq.gz
 
 
 
 
 
 
 
 
 