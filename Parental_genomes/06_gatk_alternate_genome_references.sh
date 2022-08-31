#!/bin/bash

#PBS -N gatk_alternate_reference
#PBS -l walltime=00:30:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load gatk/4.1.9.0
module load java/1.8

# Need .fai and .dict for genome before can begin
# module load samtools/1.9
# samtools faidx GCF_000214255.1_Bter_1.0_genomic.fa 
# module load picard/2.6.0
# java -jar /cm/shared/apps/picard/2.6.0/picard.jar \
# CreateSequenceDictionary R= GCF_000214255.1_Bter_1.0_genomic.fa \
# O= GCF_000214255.1_Bter_1.0_genomic.dict        

# Index the .vcf files
for file in $(ls *final.vcf)
do
    gatk IndexFeatureFile -I ${file}
done

# Run GATK to make the new alternative reference genomes
for file in $(ls *final.vcf)
do
    base=$(basename ${file} "_final.vcf")
    gatk FastaAlternateReferenceMaker \
    -R /scratch/monoallelic/hjm32/bumblebee/genome/GCF_000214255.1_Bter_1.0_genomic.fa \
    -V ${file} \
    -O ${base}.fasta 
done

# Clean up chromosome names etc
rm *fai
for file in *fasta ; do mv "${file}" "${file/.fasta/.fa}"; done

module load samtools/1.9
for file in $(ls *fa)
do
    sed -i 's/.*NC_/>NC_/g' ${file}
    sed -i 's/.*NW_/>NW_/g' ${file}
    sed -i 's/1\:1.*/1/g' ${file}
    samtools faidx ${file}
done