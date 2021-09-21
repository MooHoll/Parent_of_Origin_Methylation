#!/bin/bash

#PBS -N indel_realignment
#PBS -l walltime=08:00:00
#PBS -l vmem=40gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR

# Load software needed
module load samtools/1.9
module load picard/2.6.0 
module load gatk/3.8

echo "index bams"
for file in $(ls *.bam)
do
    java -Xms128m -Xmx1024m -Djava.io.tmpdir=$TMPDIR \
    -jar /cm/shared/apps/picard/2.6.0/picard.jar \
    BuildBamIndex \
    I=${file} \
    O=${file}.bai
done

echo "index genome and create dict file"
java -Xms128m -Xmx1024m -Djava.io.tmpdir=$TMPDIR \
    -jar /cm/shared/apps/picard/2.6.0/picard.jar CreateSequenceDictionary \
    R=honeybee_HAv3.1.fa \
    O=honeybee_HAv3.1.dict

samtools faidx honeybee_HAv3.1.fa

echo "create targets of realignment"
for file in $(ls *.bam)
do
  	base=$(basename ${file} "_sorted_RG_deduplicated.bam")
    gatk \
    -T RealignerTargetCreator \
    -R honeybee_HAv3.1.fa \
    -I ${file} \
    -o ${base}.intervals
done


echo "carry out realignment of bams"
for file in $(ls *.bam)
do
  	base=$(basename ${file} "_sorted_RG_deduplicated.bam")
    gatk \
    -T IndelRealigner \
    -R honeybee_HAv3.1.fa \
    -targetIntervals ${base}.intervals \
    -I ${file} \
    -o ${base}_indel_realigned.bam
done
