#!/bin/bash

#PBS -N alignment_lambda
#PBS -l walltime=04:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in
cd $PBS_O_WORKDIR

# Load software needed
module load bowtie2/2.3.5.1
module load samtools/1.9 

REF_FA=/scratch/monoallelic/hjm32/bumblebee/genome/lambda_genome

/scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_genome_preparation ${REF_FA}

for file in $(ls *1.fq.gz)
do
	base=$(basename $file "1.fq.gz")
	/scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark --multicore 3 --prefix lambda \
	 ${REF_FA} -1 ${base}1.fq.gz -2 ${base}2.fq.gz
done 

# Test

#!/bin/bash

#PBS -N alignment_lambda
#PBS -l walltime=04:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in
cd $PBS_O_WORKDIR

# Load software needed
module load bowtie2/2.3.5.1
module load samtools/1.9 

/scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark --multicore 3 --prefix lambda1 \
/scratch/monoallelic/hjm32/springtails/other_lambdas/lambda1 -1 AF50_WGBS_1.fq.gz -2 AF50_WGBS_2.fq.gz

/scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark --multicore 3 --prefix lambda2 \
/scratch/monoallelic/hjm32/springtails/other_lambdas/lambda2 -1 AF50_WGBS_1.fq.gz -2 AF50_WGBS_2.fq.gz

/scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark --multicore 3 --prefix lambda3 \
/scratch/monoallelic/hjm32/springtails/other_lambdas/lambda3 -1 AF50_WGBS_1.fq.gz -2 AF50_WGBS_2.fq.gz

/scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark --multicore 3 --prefix lambda4 \
/scratch/monoallelic/hjm32/springtails/other_lambdas/lambda4 -1 AF50_WGBS_1.fq.gz -2 AF50_WGBS_2.fq.gz