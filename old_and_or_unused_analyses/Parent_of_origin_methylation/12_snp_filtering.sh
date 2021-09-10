#!/bin/bash

#PBS -N filtering_snps
#PBS -l walltime=05:00:00
#PBS -l vmem=60gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load samtools/1.3.2
module load gatk/4.0.3.0 
module load java/1.8
module load picard/2.6.0

## Filter out 'fake' SNPs and poor quality calls
# By default, it filter out SNPs with quality score less than 20, 
# reads coverage more than 10, strand bias more than -0.02, 
# quality score by depth less than 1.0, mapping quality zero reads fraction more than 0.1 
# and 2 SNPs within the same 20 bp window.
for file in $(ls *.vcf)
do
  	base=$(basename ${file} "_snps.vcf")
    java -Xmx60g \
    -jar /scratch/monoallelic/hm257/bin/BisSNP-1.0.0.jar \
    -R /scratch/monoallelic/hm257/genome/GCF_000214255.1_Bter_1.0_genomic.fa \
    -T VCFpostprocess \
    -oldVcf ${file} \
    -newVcf ${base}_snp_filtered.vcf \
    -snpVcf ${file} \    
    -o ${base}_snp_filter_summary.txt
done
