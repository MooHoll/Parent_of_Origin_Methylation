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
module load gatk/3.8
module load java/1.8

# Need .fai and .dict for genome before can begin
# module load samtools/1.9
# samtools faidx GCF_000214255.1_Bter_1.0_genomic.fa 
# module load picard/2.6.0
# java -jar /cm/shared/apps/picard/2.6.0/picard.jar \
# CreateSequenceDictionary R= GCF_000214255.1_Bter_1.0_genomic.fa \
# O= GCF_000214255.1_Bter_1.0_genomic.dict

# Define file paths
REF_FA=/scratch/monoallelic/hm257/bb_imprinting/genome/GCF_000214255.1_Bter_1.0_genomic.fa                                                                                     
GATK=/cm/shared/apps/gatk/3.8/GenomeAnalysisTK.jar

# create the directory where the output files are to be written   
OUTPUT=alternate_references                                                                                                                                   
if [ ! -d "$OUTPUT" ]; then                                                                                 
    mkdir -p ${OUTPUT} 
fi                                                                                                    

# Run GATK to make the new alternative reference genomes
for file in $(ls *final.vcf)
do
    base=$(basename ${file} "_final.vcf")
    java -jar ${GATK}\
    -T FastaAlternateReferenceMaker \
    -R ${REF_FA} \
    -V ${file} \
    -o ${OUTPUT}/${base}.fasta 
done

