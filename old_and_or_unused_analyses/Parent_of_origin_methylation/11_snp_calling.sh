#############################################
# SNP calling in methylation data with BisSNP
#############################################

# Get the latest version of bis-snp onto ALICE (can re-name it as it is called 'download' but don't add file extension!)
wget https://sourceforge.net/projects/bissnp/files/latest/download
## NOTE: this can give a dodgy file, better to download onto own PC and then scp onto server

module load samtools/1.3.2
module load gatk/4.0.3.0 
module load java/1.8
module load picard/2.6.0

# Need to make sure the genome has been indexed and a .dict file created
# see: https://software.broadinstitute.org/gatk/documentation/article?id=1601
samtools faidx <genome.fa>
java -jar /cm/shared/apps/picard/2.6.0/picard.jar CreateSequenceDictionary R= <genome.fa> O= <genome.dict>

##########################################################################################

#!/bin/bash

#PBS -N sorting_bams
#PBS -l walltime=72:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load samtools/1.3.2
module load gatk/4.0.3.0 
module load java/1.8
module load picard/2.6.0

## NOTE: bam's must be sorted and have added read group information and be indexed (jeez!)
for file in $(ls *.bam)
do
  	base=$(basename ${file} "_trim_dedup_sort_RG.bam")
    java -Xmx20g \
    -jar /scratch/monoallelic/hm257/bin/BisSNP-1.0.0.jar \
    -R /scratch/monoallelic/hm257/genome/GCF_000214255.1_Bter_1.0_genomic.fa \
    -T BisulfiteGenotyper \
    -I ${file} \
    -vfn1 ${base}_cpg.vcf \
    -vfn2 ${base}_snps.vcf
done
