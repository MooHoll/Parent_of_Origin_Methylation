# Need to make one file of all the SNPs from the population
module load bcftools/1.9 
module load htslib/1.9 

# * doesn't work urgh, need a bgzip for the merge command
bgzip m08_filtered.recode.vcf
bgzip m19_filtered.recode.vcf
bgzip m23_filtered.recode.vcf
bgzip m37_filtered.recode.vcf
bgzip q08_filtered.recode.vcf
bgzip q19_filtered.recode.vcf
bgzip q23_filtered.recode.vcf
bgzip q37_filtered.recode.vcf

# Need index for the merge command
bcftools index m08_filtered.recode.vcf.gz
bcftools index m19_filtered.recode.vcf.gz
bcftools index m23_filtered.recode.vcf.gz
bcftools index m37_filtered.recode.vcf.gz
bcftools index q08_filtered.recode.vcf.gz
bcftools index q19_filtered.recode.vcf.gz
bcftools index q23_filtered.recode.vcf.gz
bcftools index q37_filtered.recode.vcf.gz

#!/bin/bash

#PBS -N merge_snps
#PBS -l walltime=00:10:00
#PBS -l vmem=100gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=3

# Run script in the working directory it was submitted in  
cd $PBS_O_WORKDIR 

module load bcftools/1.9 
module load htslib/1.9 

bcftools merge --threads 3 -o bumblebee_pop_snps.vcf \
m08_filtered.recode.vcf.gz m19_filtered.recode.vcf.gz m23_filtered.recode.vcf.gz m37_filtered.recode.vcf.gz \
q08_filtered.recode.vcf.gz q19_filtered.recode.vcf.gz q23_filtered.recode.vcf.gz q37_filtered.recode.vcf.gz

#-------------------------------------------------------

# Running SNPeff to get syn and nonsyn subsitution info
module load snpeff/4.3t 
cp /cm/shared/apps/snpeff/4.3t/snpEff.config.example ./snpEff.config

# Need to add the bumblebee genome information to the config file
nano snpEff.config
# Added under non-standard databases:

# Bumblebee genome: version Bter_1.0
Bter_1.0.genome : Bumblebee

# now make the database
mkdir data
cd data
mkdir Bter_1.0
mkdir genomes

mv ref_Bter_1.0_top_level.gff3 genes.gff
mv genes.gff data/Bter_1.0/

mv GCF_000214255.1_Bter_1.0_genomic.fa Bter_1.0.fa
mv Bter_1.0.fa data/genomes/

snpEff build -c snpEff.config -gff3 -v Bter_1.0

# Following andrew's powerpoint, filter out intergenic varients and deal with alt transcripts
snpEff -c snpEff.config -v -no-downstream -no-intron -no-upstream -no-utr -no-intergenic \
-stats bumblebee_pop_filterstats.html -canon Bter_1.0 \
bumblebee_pop_snps.vcf > bumblebee_pop_filtered_for_pnps.vcf

#SnpSift isn't avaliable through the ALICE module for some reason so download manually

# Download latest version
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
# Unzip file
unzip snpEff_latest_core.zip

# Extract information
module load java/11.0.2

java -jar /scratch/monoallelic/hjm32/bin/snpeff/snpEff/SnpSift.jar extractFields \
bumblebee_pop_filtered_for_pnps.vcf CHROM POS REF ALT ANN[*].EFFECT > bumblebee_pop_snps_with_info.txt

# Now need to count up synonymous and missense (nonsyn) SNP annotations per gene, after labelling each
# SNP with it's gene_id