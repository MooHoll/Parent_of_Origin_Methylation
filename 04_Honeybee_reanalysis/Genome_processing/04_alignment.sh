#!/bin/bash

#PBS -N alignment
#PBS -l walltime=00:05:00
#PBS -l vmem=40gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16
#PBS -q devel

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR

# Load software needed
module load bowtie2/2.3.5.1
module load samtools/1.9

# Prep genome
bowtie2-build /scratch/monoallelic/hjm31/honeybee_reanalsysis/genomes/honeybee/honeybee_HAv3.1.fa honeybee

BOWTIE2_INDEXES=/scratch/monoallelic/hjm31/honeybee_reanalsysis/genomes/honeybee

echo "alignment"
for file in $(ls *_1.fq.gz)
do
	base=$(basename $file "_1.fq.gz")
    bowtie2 --sensitive --threads 14 \
    -x honeybee \
    -1 ${base}_1.fq.gz \
    -2 ${base}_2.fq.gz \
    -U ${base}_unpaired_1.fq.gz,${base}_unpaired_2.fq.gz \
    -S ${base}.sam
done


echo "convert sams to bams"
for file in $(ls *.sam)
do
	base=$(basename $file ".sam")
    samtools view -bS ${base}.sam > ${base}.bam
done

echo "sort bams"
for file in $(ls *.bam)
do
	base=$(basename $file ".bam")
    samtools sort -@ 14 -o ${base}_sorted.bam ${file}
done