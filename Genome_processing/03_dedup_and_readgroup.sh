#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/bombus_terrestris/logs/dedup_and_readgroup_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/misc/analyses/bombus_terrestris/bams/*sorted.bam ./

#---------------------------------------------

echo "indexing bams"
for file in $(ls *.bam)
do
samtools index ${file}
done

echo "adding read groups"
for file in $(ls *sorted.bam)
do
  	base=$(basename ${file} "_sorted.bam")
    picard \
    AddOrReplaceReadGroups \
    I=${file} \
    O=${base}_sorted_RG.bam \
    RGID=0001${base} \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=NA \
    RGSM=${base}
done

echo "adding read groups"
for file in $(ls *RG.bam)
do
  	base=$(basename ${file} "_sorted_RG.bam")
    picard \
    MarkDuplicates \
    I=${file} \
    O=${base}_sorted_RG_deduplicated.bam \
    REMOVE_DUPLICATES=TRUE
done

#---------------------------------------------

echo "moving outputs"
mv *_sorted_RG_deduplicated.bam /data/ross/misc/analyses/bombus_terrestris/bams

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"