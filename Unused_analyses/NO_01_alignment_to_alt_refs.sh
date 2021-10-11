# Make a folder for each new reference genome and place .fa file inside

# Also need to edit the new genome fastas to make the chromosome names match all other files generated
sed -i 's/>.*NC/>NC/g' m08.fasta
sed -i 's/>.*NW/>NW/g' m08.fasta
sed -i 's/:1//g' m08.fasta

sed -i 's/>.*NC/>NC/g' m19.fasta
sed -i 's/>.*NW/>NW/g' m19.fasta
sed -i 's/:1//g' m19.fasta

sed -i 's/>.*NC/>NC/g' m23.fasta
sed -i 's/>.*NW/>NW/g' m23.fasta
sed -i 's/:1//g' m23.fasta

sed -i 's/>.*NC/>NC/g' m37.fasta
sed -i 's/>.*NW/>NW/g' m37.fasta
sed -i 's/:1//g' m37.fasta

sed -i 's/>.*NC/>NC/g' q08.fasta
sed -i 's/>.*NW/>NW/g' q08.fasta
sed -i 's/:1//g' q08.fasta

sed -i 's/>.*NC/>NC/g' q19.fasta
sed -i 's/>.*NW/>NW/g' q19.fasta
sed -i 's/:1//g' q19.fasta

sed -i 's/>.*NC/>NC/g' q23.fasta
sed -i 's/>.*NW/>NW/g' q23.fasta
sed -i 's/:1//g' q23.fasta

sed -i 's/>.*NC/>NC/g' q37.fasta
sed -i 's/>.*NW/>NW/g' q37.fasta
sed -i 's/:1//g' q37.fasta

module load bowtie2/2.3.5.1 
module load samtools/1.9

for folder in $(ls -d *_genome)
do
	/scratch/monoallelic/hjm31/bin_new/Bismark-0.22.3/bismark_genome_preparation ./${folder}
done

# -----------------------------------
# Alignment to parental genomes
# -----------------------------------

#!/bin/bash

#PBS -N alignment_to_ref
#PBS -l walltime=25:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=20

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load bowtie2/2.3.5.1 
module load samtools/1.9


/scratch/monoallelic/hjm31/bin_new/Bismark-0.22.3/bismark --multicore 3 --prefix male --non_bs_mm \
/scratch/monoallelic/hjm31/bumblebee/genome/alternate_references/m08_genome \
 -1 trim_w08_1.fq.gz -2 trim_w08_2.fq.gz
 
/scratch/monoallelic/hjm31/bin_new/Bismark-0.22.3/bismark --multicore 3 --prefix queen --non_bs_mm \
/scratch/monoallelic/hjm31/bumblebee/genome/alternate_references/q08_genome \
 -1 trim_w08_1.fq.gz -2 trim_w08_2.fq.gz
 

/scratch/monoallelic/hjm31/bin_new/Bismark-0.22.3/bismark --multicore 3 --prefix male --non_bs_mm \
/scratch/monoallelic/hjm31/bumblebee/genome/alternate_references/m19_genome \
 -1 trim_w19_1.fq.gz -2 trim_w19_2.fq.gz
 
/scratch/monoallelic/hjm31/bin_new/Bismark-0.22.3/bismark --multicore 3 --prefix queen --non_bs_mm \
/scratch/monoallelic/hjm31/bumblebee/genome/alternate_references/q19_genome \
 -1 trim_w19_1.fq.gz -2 trim_w19_2.fq.gz
 

/scratch/monoallelic/hjm31/bin_new/Bismark-0.22.3/bismark --multicore 3 --prefix male --non_bs_mm \
/scratch/monoallelic/hjm31/bumblebee/genome/alternate_references/m23_genome \
 -1 trim_w23_1.fq.gz -2 trim_w23_2.fq.gz
 
/scratch/monoallelic/hjm31/bin_new/Bismark-0.22.3/bismark --multicore 3 --prefix queen --non_bs_mm \
/scratch/monoallelic/hjm31/bumblebee/genome/alternate_references/q23_genome \
 -1 trim_w23_1.fq.gz -2 trim_w23_2.fq.gz

 
/scratch/monoallelic/hjm31/bin_new/Bismark-0.22.3/bismark --multicore 3 --prefix male --non_bs_mm \
/scratch/monoallelic/hjm31/bumblebee/genome/alternate_references/m37_genome \
 -1 trim_w37_1.fq.gz -2 trim_w37_2.fq.gz
 
/scratch/monoallelic/hjm31/bin_new/Bismark-0.22.3/bismark --multicore 3 --prefix queen --non_bs_mm \
/scratch/monoallelic/hjm31/bumblebee/genome/alternate_references/q37_genome \
 -1 trim_w37_1.fq.gz -2 trim_w37_2.fq.gz
 
 
for file in $(ls *.bam)
do
    /scratch/monoallelic/hjm31/bin_new/Bismark-0.22.3/deduplicate_bismark -p ${file}
done
 
 
 
 
 
 