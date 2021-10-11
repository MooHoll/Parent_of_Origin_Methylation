# -----------------------------------------------
# Get queen homozygous alternative SNPS: honeybee
# -----------------------------------------------

# Count total heterozygote SNPs (queens only)
grep -c 0/1 <sample.vcf>
# 875 537001
# 882 602919
# 888 1242218
# 894 937190

# Count total homozygous SNPs (queens only)
grep -c 1/1 <sample.vcf>
# 875 378275
# 882 396595
# 888 672170
# 894 488673

# Write file containing only the alternate homozygous SNPs (queens only)
grep -e ^# -e "1/1" 875_queen_filtered.recode.vcf > 875_queen_alt_homozygous_snps.vcf
grep -e ^# -e "1/1" 882_queen_filtered.recode.vcf > 882_queen_alt_homozygous_snps.vcf
grep -e ^# -e "1/1" 888_queen_filtered.recode.vcf > 888_queen_alt_homozygous_snps.vcf
grep -e ^# -e "1/1" 894_queen_filtered.recode.vcf > 894_queen_alt_homozygous_snps.vcf

# Count the male SNPs (all already alternative homozygous as haploid)
grep -v ^# <sample.vcf> | wc -l
# 875 1200722
# 882 1236225
# 888 738195
# 894 951432

# -----------------------------------------------
# Make files with SNPs unique to each parent
# -----------------------------------------------

module load bedtools/2.28.0

# Remove uninformative SNPs that males and queens have in common
# -A needed to remove whole SNP information if there is overlap, -header needed for future commands to recognise the .vcf format 
subtractBed -header -A -a 875_drone_filtered.recode.vcf -b 875_queen_alt_homozygous_snps.vcf > 875_drone_unique.vcf
grep -v ^# 875_drone_unique.vcf | wc -l
# 985500

subtractBed -header -A -a 875_queen_alt_homozygous_snps.vcf -b 875_drone_filtered.recode.vcf > 875_queen_unique.vcf
grep -v ^# 875_queen_unique.vcf | wc -l
# 163053

subtractBed -header -A -a 882_drone_filtered.recode.vcf -b 882_queen_alt_homozygous_snps.vcf > 882_drone_unique.vcf
grep -v ^# 882_drone_unique.vcf | wc -l
# 1007359

subtractBed -header -A -a 882_queen_alt_homozygous_snps.vcf -b 882_drone_filtered.recode.vcf > 882_queen_unique.vcf
grep -v ^# 882_queen_unique.vcf | wc -l
# 167729

subtractBed -header -A -a 888_drone_filtered.recode.vcf -b 888_queen_alt_homozygous_snps.vcf > 888_drone_unique.vcf
grep -v ^# 888_drone_unique.vcf | wc -l
# 506936

subtractBed -header -A -a 888_queen_alt_homozygous_snps.vcf -b 888_drone_filtered.recode.vcf > 888_queen_unique.vcf
grep -v ^# 888_queen_unique.vcf | wc -l
# 440911

subtractBed -header -A -a 894_drone_filtered.recode.vcf -b 894_queen_alt_homozygous_snps.vcf > 894_drone_unique.vcf
grep -v ^# 894_drone_unique.vcf | wc -l
# 716899

subtractBed -header -A -a 894_queen_alt_homozygous_snps.vcf -b 894_drone_filtered.recode.vcf > 894_queen_unique.vcf
grep -v ^# 894_queen_unique.vcf | wc -l
# 254140

# Now filter out all the CT and TC SNPs
awk '!($4 == "C" && $5 == "T") {print ;}' 875_queen_unique.vcf > inbetween
awk '!($4 == "T" && $5 == "C") {print ;}' inbetween > 875_queen_final.vcf
grep -v ^# 875_queen_final.vcf | wc -l
# 97509

awk '!($4 == "C" && $5 == "T") {print ;}' 882_queen_unique.vcf > inbetween
awk '!($4 == "T" && $5 == "C") {print ;}' inbetween > 882_queen_final.vcf
grep -v ^# 882_queen_final.vcf | wc -l
# 99117

awk '!($4 == "C" && $5 == "T") {print ;}' 888_queen_unique.vcf > inbetween
awk '!($4 == "T" && $5 == "C") {print ;}' inbetween > 888_queen_final.vcf
grep -v ^# 888_queen_final.vcf | wc -l
# 263615

awk '!($4 == "C" && $5 == "T") {print ;}' 894_queen_unique.vcf > inbetween
awk '!($4 == "T" && $5 == "C") {print ;}' inbetween > 894_queen_final.vcf
grep -v ^# 894_queen_final.vcf | wc -l
# 151845


awk '!($4 == "C" && $5 == "T") {print ;}' 875_drone_unique.vcf > inbetween
awk '!($4 == "T" && $5 == "C") {print ;}' inbetween > 875_drone_final.vcf
grep -v ^# 875_drone_final.vcf | wc -l
# 593515

awk '!($4 == "C" && $5 == "T") {print ;}' 882_drone_unique.vcf > inbetween
awk '!($4 == "T" && $5 == "C") {print ;}' inbetween > 882_drone_final.vcf
grep -v ^# 882_drone_final.vcf | wc -l
# 606131

awk '!($4 == "C" && $5 == "T") {print ;}' 888_drone_unique.vcf > inbetween
awk '!($4 == "T" && $5 == "C") {print ;}' inbetween > 888_drone_final.vcf
grep -v ^# 888_drone_final.vcf | wc -l
# 304708

awk '!($4 == "C" && $5 == "T") {print ;}' 894_drone_unique.vcf > inbetween
awk '!($4 == "T" && $5 == "C") {print ;}' inbetween > 894_drone_final.vcf
grep -v ^# 894_drone_final.vcf | wc -l
# 432261