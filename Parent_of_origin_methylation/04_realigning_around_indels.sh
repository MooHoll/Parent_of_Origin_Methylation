#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/bombus_terrestris/logs/realignment_around_indels_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/misc/analyses/bombus_terrestris/bams/*_sorted_RG_deduplicated.bam ./
rsync /data/ross/misc/analyses/bombus_terrestris/genome/GCF_000214255.1_Bter_1.0_genomic.fa ./
#---------------------------------------------

echo "index bams"
for file in $(ls *.bam)
do
    picard \
    BuildBamIndex \
    I=${file} \
    O=${file}.bai
done

echo "index genome and create dict file"
picard CreateSequenceDictionary -R GCF_000214255.1_Bter_1.0_genomic.fa
samtools faidx GCF_000214255.1_Bter_1.0_genomic.fa

echo "create targets of realignment"
for file in $(ls *.bam)
do
  	base=$(basename ${file} "_sorted_RG_deduplicated.bam")
    gatk3 \
    -T RealignerTargetCreator \
    -R GCF_000214255.1_Bter_1.0_genomic.fa \
    -I ${file} \
    -o ${base}.intervals
done


echo "carry out realignment of bams"
for file in $(ls *.bam)
do
  	base=$(basename ${file} "_sorted_RG_deduplicated.bam")
    gatk3 \
    -T IndelRealigner \
    -R GCF_000214255.1_Bter_1.0_genomic.fa \
    -targetIntervals ${base}.intervals \
    -I ${file} \
    -o ${base}_indel_realigned.bam
done


#---------------------------------------------

echo "moving outputs"
mv *_indel_realigned.bam /data/ross/misc/analyses/bombus_terrestris/bams

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"