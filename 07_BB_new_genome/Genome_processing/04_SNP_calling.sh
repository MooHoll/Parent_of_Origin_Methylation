#!/bin/bash

#PBS -N calling_SNPs
#PBS -l walltime=00:05:00
#PBS -l vmem=40gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8
#PBS -q devel

# Run script in the working directory it was submitted in (8hrs)
cd $PBS_O_WORKDIR

# Load software needed
module load samtools/1.9
module load freebayes/1.1.0
module load vcftools/0.1.14

males=m*.bam
queens=q*.bam

echo "call SNPs on queens"
# min count 2 reads for alternative alleles, 
# min 5 reads per SNP, ignore complex events, indels and mnps
for file in ${queens}
do
  	base=$(basename ${file} "_indel_realigned.bam")
    freebayes \
        -f GCA_910591885.1.fa \
        -C 2 \
        --min-coverage 5 \
        -u \
        -i \
        -X \
        -b ${file} \
        > ${base}_vcf.gz
done

echo "call SNPs on males"
# ploidy set to 1
for file in ${males}
do
    base=$(basename ${file} "_indel_realigned.bam")
    freebayes \
        -f GCA_910591885.1.fa \
        -p 1 \
        --min-coverage 5 \
        -C 2 \
        -u \
        -i \
        -X \
        -b ${file} \
        > ${base}_vcf.gz
done

for file in $(ls *gz)
do
	base=$(basename ${file} "vcf.gz")
    vcftools --gzvcf ${file} --max-alleles 2 --minQ 20 --min-meanDP 10 \
    --recode --recode-INFO-all --out ${base}filtered
done
