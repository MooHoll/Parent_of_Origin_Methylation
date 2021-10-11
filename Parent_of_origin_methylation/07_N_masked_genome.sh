# Make an N masked genome for each colony

module load bedtools/2.28.0
module load bcftools/1.9
module load htslib/1.9 
module load samtools/1.9

# First need one SNP file for each colony which contains both the unique mother
# and father SNPs so they all get masked

# * doesn't work urgh, need a bgzip for the merge command
bgzip 875_drone_final.vcf
bgzip 882_drone_final.vcf
bgzip 888_drone_final.vcf
bgzip 894_drone_final.vcf
bgzip 875_queen_final.vcf
bgzip 882_queen_final.vcf
bgzip 888_queen_final.vcf
bgzip 894_queen_final.vcf

# Need index for the merge command
bcftools index 875_drone_final.vcf.gz
bcftools index 882_drone_final.vcf.gz
bcftools index 888_drone_final.vcf.gz
bcftools index 894_drone_final.vcf.gz
bcftools index 875_queen_final.vcf.gz
bcftools index 882_queen_final.vcf.gz
bcftools index 888_queen_final.vcf.gz
bcftools index 894_queen_final.vcf.gz

# Make one SNP file per colony
bcftools merge -o 875_colony.vcf.gz 875_drone_final.vcf.gz 875_queen_final.vcf.gz
bcftools merge -o 882_colony.vcf.gz 882_drone_final.vcf.gz 882_queen_final.vcf.gz
bcftools merge -o 888_colony.vcf.gz 888_drone_final.vcf.gz 888_queen_final.vcf.gz
bcftools merge -o 894_colony.vcf.gz 894_drone_final.vcf.gz 894_queen_final.vcf.gz

# Make an N masked genome per colony
bedtools maskfasta -fi honeybee_HAv3.1.fa -bed ../../snps/875_colony.vcf.gz -fo 875_masked_genome.fa 
bedtools maskfasta -fi honeybee_HAv3.1.fa -bed ../../snps/882_colony.vcf.gz -fo 882_masked_genome.fa 
bedtools maskfasta -fi honeybee_HAv3.1.fa -bed ../../snps/888_colony.vcf.gz -fo 888_masked_genome.fa 
bedtools maskfasta -fi honeybee_HAv3.1.fa -bed ../../snps/894_colony.vcf.gz -fo 894_masked_genome.fa 