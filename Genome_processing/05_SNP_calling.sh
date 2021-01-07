#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/bombus_terrestris/logs/SNPcalling_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/misc/analyses/bombus_terrestris/bams/*_indel_realigned.bam ./
rsync /data/ross/misc/analyses/bombus_terrestris/genome/GCF_000214255.1_Bter_1.0_genomic.fa ./

#---------------------------------------------

# NOTE: freebayes is single thread, can multithread with freebayes-parallel
# but this only splits the genome using a simple python script, will just avoid

echo "indexing bams"
for file in $(ls *.bam)
do
samtools index ${file}
done

males=m*.bam
queens=q*.bam

echo "call SNPs on queens"
# min count 2 reads for alternative alleles, 
# min 5 reads per SNP, ignore complex events, indels and mnps
for file in ${queens}
do
  	base=$(basename ${file} "_indel_realigned.bam")
    freebayes \
        -f GCF_000214255.1_Bter_1.0_genomic.fa \
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
        -f GCF_000214255.1_Bter_1.0_genomic.fa \
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
    vcftools --gzvcf ${file} --max-alleles 2 --minQ 20 --min-meanDP 10 --recode --recode-INFO-all --out ${base}filtered
done

#---------------------------------------------

echo "moving outputs"
mv * /data/ross/misc/analyses/bombus_terrestris/snps

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"