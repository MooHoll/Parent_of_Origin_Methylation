#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/bombus_terrestris/logs/fastqc_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in from softlinks"
rsync -L /data/ross/misc/analyses/bombus_terrestris/raw/*.fq.gz ./

echo "doing the shiz"
for file in $(ls *.gz)
do
	fastqc -t 9 ${file}
done

echo "moving outputs"
mv *html /data/ross/misc/analyses/bombus_terrestris/raw

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"