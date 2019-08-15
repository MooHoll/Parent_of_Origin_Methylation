### Ran by hand on the command line to count SNPs and remove non-informative queen SNPs

# Count total SNPs in a file
grep -v '#' <file.vcf> | wc -l

# Count total heterozygote SNPs (queens only)
zgrep -v '#' <sample.vcf.gz> | zgrep -vc "1[/|]1"

# Count total homozygous SNPs (queens only)
zgrep -v '#' <sample.vcf.gz> | zgrep -c "1[/|]1" 

# Write file containing only the alternate homozygous SNPs (queens only)
# Do the double negative to keep headers etc!
zgrep -v "0[/|]1" <sample.vcf.gz> | zgrep -v "1[/|]2" > <output.vcf.gz>

## There are no 0/0 dufus as it means there is no SNP >.<
## But there are 1/2 which means heterozygous but two alternate SNPs

#### NOTE: as the SNP caller used does account for ploidy also need to filter male SNPs to keep homo-alt
### some mistakes in SNP calling mean we have some 'heterozygous' points.



# Ran by hand on the commandline to produce .vcf files with SNPs unique to either the mother or father 

module load bedtools/2.25.0

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














