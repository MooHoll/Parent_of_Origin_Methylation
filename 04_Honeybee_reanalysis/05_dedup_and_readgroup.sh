#!/bin/bash

#PBS -N read_groups
#PBS -l walltime=00:05:00
#PBS -l vmem=40gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8
#PBS -q devel

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR

# Load software needed
module load samtools/1.9
module load picard/2.6.0 

#echo "indexing bams"
#for file in $(ls *sorted.bam)
#do
#  samtools index ${file}
#done

echo "adding read groups"
for file in $(ls *sorted.bam)
do
  	base=$(basename ${file} "_trim_sorted.bam")
    java -Xms128m -Xmx1024m -Djava.io.tmpdir=$TMPDIR \
    -jar /cm/shared/apps/picard/2.6.0/picard.jar \
    AddOrReplaceReadGroups \
    I=${file} \
    O=${base}_sorted_RG.bam \
    RGID=0001${base} \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=NA \
    RGSM=${base}
done

echo "deduplicating"
for file in $(ls *RG.bam)
do
  	base=$(basename ${file} "_sorted_RG.bam")
    java -Xms128m -Xmx1024m -Djava.io.tmpdir=$TMPDIR \
    -jar /cm/shared/apps/picard/2.6.0/picard.jar \
    MarkDuplicates \
    I=${file} \
    O=${base}_sorted_RG_deduplicated.bam \
    M=${base}_duplicate_metrics.txt \
    REMOVE_DUPLICATES=TRUE
done