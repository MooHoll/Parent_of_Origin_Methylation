#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/bombus_terrestris/logs/alignment_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------
# Alignment to reference genome

echo "copying data in"
rsync -L /data/ross/misc/analyses/bombus_terrestris/raw/*.fq.gz ./
rsync /data/ross/misc/analyses/bombus_terrestris/genome/GCF_000214255.1_Bter_1.0_genomic.fa ./

echo "making the genome index"
bowtie2-build ./GCF_000214255.1_Bter_1.0_genomic.fa b_terrestris

echo "starting alignment"
for file in $(ls *_1.fq.gz)
do
	base=$(basename $file "_1.fq.gz")
    bowtie2 --sensitive --threads 30 \
    -x b_terrestris \
    -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz \
    -S ${base}.sam
done

# consider adding --trim5 with 5bp

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
    samtools sort -@ 25 -o ${base}_sorted.bam ${file}
done

echo "moving outputs"
mv *sorted.bam /data/ross/misc/analyses/bombus_terrestris

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"