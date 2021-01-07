# Ran by hand on the commandline to produce .vcf files with SNPs unique to either the mother or father 

module load bedtools/2.25.0

# Rename files
for file in m*.vcf; do mv "${file}" "${file/_filtered.recode/}"; done
for file in q*.vcf; do mv "${file}" "${file/_filtered_final/}"; done

# Remove uninformative SNPs that males and queens have in common
# -A needed to remove whole SNP information if there is overlap, -header needed for future commands to recognise the .vcf format 
subtractBed -header -A -a m02.vcf -b q02.vcf > m02_unique.vcf
subtractBed -header -A -a q02.vcf -b m02.vcf > q02_unique.vcf

subtractBed -header -A -a m12.vcf -b q12.vcf > m12_unique.vcf
subtractBed -header -A -a q12.vcf -b m12.vcf > q12_unique.vcf

subtractBed -header -A -a m22.vcf -b q22.vcf > m22_unique.vcf
subtractBed -header -A -a q22.vcf -b m22.vcf > q22_unique.vcf

subtractBed -header -A -a m31.vcf -b q31.vcf > m31_unique.vcf
subtractBed -header -A -a q31.vcf -b m31.vcf > q31_unique.vcf

# Make a file of common SNPs for all crosses
intersectBed -header -a m02_unique.vcf -b m12_unique.vcf m22_unique.vcf m31_unique.vcf > male_SNPs.vcf
intersectBed -header -a q02_unique.vcf -b q12_unique.vcf q22_unique.vcf q31_unique.vcf > queen_SNPs.vcf
