# Make an N masked genome for each colony

module load bedtools/2.28.0
module load bcftools/1.9
module load htslib/1.9 
module load samtools/1.9

# First need one SNP file for each colony which contains both the unique mother
# and father SNPs so they all get masked

# * doesn't work urgh, need a bgzip for the merge command
bgzip m08_final.vcf
bgzip m19_final.vcf
bgzip m23_final.vcf
bgzip m37_final.vcf
bgzip q08_final.vcf
bgzip q19_final.vcf
bgzip q23_final.vcf
bgzip q37_final.vcf

# Need index for the merge command
bcftools index *.vcf.gz

# Make one SNP file per colony
bcftools merge -o 08_colony.vcf.gz m08_final.vcf.gz q08_final.vcf.gz
bcftools merge -o 19_colony.vcf.gz m19_final.vcf.gz q19_final.vcf.gz
bcftools merge -o 23_colony.vcf.gz m23_final.vcf.gz q23_final.vcf.gz
bcftools merge -o 37_colony.vcf.gz m37_final.vcf.gz q37_final.vcf.gz

# Make an N masked genome per colony
bedtools maskfasta -fi GCF_000214255.1_Bter_1.0_genomic.fa -bed ../snps/08_colony.vcf.gz -fo 08_masked_genome.fa 
bedtools maskfasta -fi GCF_000214255.1_Bter_1.0_genomic.fa -bed ../snps/19_colony.vcf.gz -fo 19_masked_genome.fa 
bedtools maskfasta -fi GCF_000214255.1_Bter_1.0_genomic.fa -bed ../snps/23_colony.vcf.gz -fo 23_masked_genome.fa 
bedtools maskfasta -fi GCF_000214255.1_Bter_1.0_genomic.fa -bed ../snps/37_colony.vcf.gz -fo 37_masked_genome.fa 