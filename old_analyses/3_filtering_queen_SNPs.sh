# Ran by hand on the command line to count SNPs and remove non-informative queen SNPs

# Keep quality over 20
vcftools --gzvcf <vcf> --out <sample>filtered --minQ 20 --recode --recode-INFO-all

# Count total heterozygote SNPs (queens only)
zgrep -v ^# <sample.vcf.gz> | zgrep -c "0[/|]1"

# Count total homozygous SNPs (queens only)
zgrep -v ^# <sample.vcf.gz> | zgrep -vc "0[/|]1" 

# Write file containing only the alternate homozygous SNPs (queens only)
zgrep -v "0[/|]1" <sample.vcf.gz> | zgrep -v "1[/|]2" > <output.vcf.gz>
