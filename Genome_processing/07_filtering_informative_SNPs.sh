# -----------------------------------------------
# Get queen homozygous alternative SNPS
# -----------------------------------------------

# Count total heterozygote SNPs (queens only)
grep -c 0/1 <sample.vcf>
# q08 608432
# q19 556697
# q23 578448
# q37 561185

# Count total homozygous SNPs (queens only)
grep -c 1/1 <sample.vcf>
# q08 409964
# q19 379222
# q23 350774
# q37 362659

# Write file containing only the alternate homozygous SNPs (queens only)
grep -e ^# -e "1/1" q08_filtered.recode.vcf > q08_alt_homozygous_snps.vcf
grep -e ^# -e "1/1" q19_filtered.recode.vcf > q19_alt_homozygous_snps.vcf
grep -e ^# -e "1/1" q23_filtered.recode.vcf > q23_alt_homozygous_snps.vcf
grep -e ^# -e "1/1" q37_filtered.recode.vcf > q37_alt_homozygous_snps.vcf

# Count the male SNPs (all already alternative homozygous as haploid)
grep -v ^# <sample.vcf> | wc -l
# m08 732610
# m19 676705
# m23 663441
# m37 743465

# -----------------------------------------------
# Make files with SNPs unique to each parent
# -----------------------------------------------

module load bedtools/2.28.0

# Remove uninformative SNPs that males and queens have in common
# -A needed to remove whole SNP information if there is overlap, -header needed for future commands to recognise the .vcf format 
subtractBed -header -A -a m08_filtered.recode.vcf -b q08_alt_homozygous_snps.vcf > m08_unique.vcf
grep -v ^# m08_unique.vcf | wc -l
# 456283

subtractBed -header -A -a q08_alt_homozygous_snps.vcf -b m08_filtered.recode.vcf > q08_unique.vcf
grep -v ^# q08_unique.vcf | wc -l
# 133637

subtractBed -header -A -a m19_filtered.recode.vcf -b q19_alt_homozygous_snps.vcf > m19_unique.vcf
grep -v ^# m19_unique.vcf | wc -l
# 434965

subtractBed -header -A -a q19_alt_homozygous_snps.vcf -b m19_filtered.recode.vcf > q19_unique.vcf
grep -v ^# q19_unique.vcf | wc -l
# 137482

subtractBed -header -A -a m23_filtered.recode.vcf -b q23_alt_homozygous_snps.vcf > m23_unique.vcf
grep -v ^# m23_unique.vcf | wc -l
# 443122

subtractBed -header -A -a q23_alt_homozygous_snps.vcf -b m23_filtered.recode.vcf > q23_unique.vcf
grep -v ^# q23_unique.vcf | wc -l
# 130455

subtractBed -header -A -a m37_filtered.recode.vcf -b q37_alt_homozygous_snps.vcf > m37_unique.vcf
grep -v ^# m37_unique.vcf | wc -l
# 498411

subtractBed -header -A -a q37_alt_homozygous_snps.vcf -b m37_filtered.recode.vcf > q37_unique.vcf
grep -v ^# q37_unique.vcf | wc -l
# 117605