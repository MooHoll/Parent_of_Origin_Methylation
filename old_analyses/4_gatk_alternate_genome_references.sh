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

# Define file paths
REF_FA=/scratch/monoallelic/hm257/bb_imprinting/genome                                                                                     
GATK=/cm/shared/apps/gatk/3.8/GenomeAnalysisTK.jar

# create the directory where the output files are to be written   
OUTPUT=alternate_references                                                                                                                                   
if [ ! -d "$OUTPUT" ]; then                                                                                 
    mkdir -p ${OUTPUT} 
fi


# Run GATK to make the new alternative reference genomes
for file in $(ls *.vcf)
do
    base=$(basename ${file} "_filtered.recode.vcf")
    java -jar ${GATK}\
    -T FastaAlternateReferenceMaker \
    -R ${REF_FA} \
    -V ${file} \
    -o ${OUTPUT}/${base}.fasta 
done

