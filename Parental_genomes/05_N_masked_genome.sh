# Make an N masked genome for each colony

module load bedtools/2.28.0
module load bcftools/1.9
module load htslib/1.9 
module load samtools/1.9

# First need one SNP file for each colony which contains both the mother
# and father SNPs so they all get masked

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

# Make one SNP file per colony
bcftools merge -o 08_colony.vcf.gz m08_filtered.recode.vcf.gz q08_filtered.recode.vcf.gz
bcftools merge -o 19_colony.vcf.gz m19_filtered.recode.vcf.gz q19_filtered.recode.vcf.gz
bcftools merge -o 23_colony.vcf.gz m23_filtered.recode.vcf.gz q23_filtered.recode.vcf.gz
bcftools merge -o 37_colony.vcf.gz m37_filtered.recode.vcf.gz q37_filtered.recode.vcf.gz

# Make an N masked genome per colony
bedtools maskfasta -fi GCF_000214255.1_Bter_1.0_genomic.fa -bed 08_colony.vcf.gz -fo 08_masked_genome.fa 
bedtools maskfasta -fi GCF_000214255.1_Bter_1.0_genomic.fa -bed 19_colony.vcf.gz -fo 19_masked_genome.fa 
bedtools maskfasta -fi GCF_000214255.1_Bter_1.0_genomic.fa -bed 23_colony.vcf.gz -fo 23_masked_genome.fa 
bedtools maskfasta -fi GCF_000214255.1_Bter_1.0_genomic.fa -bed 37_colony.vcf.gz -fo 37_masked_genome.fa 