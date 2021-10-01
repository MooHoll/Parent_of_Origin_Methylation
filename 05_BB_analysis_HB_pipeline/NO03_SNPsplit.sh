# Make SNP files for SNP split to work
module load htslib/1.9
bgzip -d *_final.vcf.gz

for file in $(ls *_final.vcf)
do
    base=$(basename ${file} "_final.vcf")
    grep -v '#' ${base}_final.vcf |
    awk '{print $6 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "/" $5}' - > ${base}_1.txt
    echo -e "ID\tChr\tPosition\tSNP value\tRef/SNP" | cat - ${base}_1.txt > ${base}_snpsplit.txt
done

rm *1.txt

#---------------------------------------------------------

#!/bin/bash

#PBS -N SNP_split
#PBS -l walltime=08:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=2

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load samtools/1.9

# Run snpsplit for each parental SNP set
/scratch/monoallelic/hjm31/bin_new/SNPsplit-0.5.0/SNPsplit \
--paired --bisulfite -o ./male_08 \
--snp_file ../snps/m08_snpsplit.txt \
trim_w08_1_bismark_bt2_pe.deduplicated.bam
 
 /scratch/monoallelic/hjm31/bin_new/SNPsplit-0.5.0/SNPsplit \
--paired --bisulfite -o ./queen_08 \
--snp_file ../snps/q08_snpsplit.txt \
trim_w08_1_bismark_bt2_pe.deduplicated.bam

#---

/scratch/monoallelic/hjm31/bin_new/SNPsplit-0.5.0/SNPsplit \
--paired --bisulfite ./male_19 \
--snp_file ../snps/m19_snpsplit.txt \
trim_w19_1_bismark_bt2_pe.deduplicated.bam
 
 /scratch/monoallelic/hjm31/bin_new/SNPsplit-0.5.0/SNPsplit \
--paired --bisulfite -o ./queen_19 \
--snp_file ../snps/q19_snpsplit.txt \
trim_w19_1_bismark_bt2_pe.deduplicated.bam

#---

/scratch/monoallelic/hjm31/bin_new/SNPsplit-0.5.0/SNPsplit \
--paired --bisulfite ./male_23 \
--snp_file ../snps/m23_snpsplit.txt \
trim_w23_1_bismark_bt2_pe.deduplicated.bam
 
 /scratch/monoallelic/hjm31/bin_new/SNPsplit-0.5.0/SNPsplit \
--paired --bisulfite -o ./queen_23 \
--snp_file ../snps/q23_snpsplit.txt \
trim_w23_1_bismark_bt2_pe.deduplicated.bam

#---

/scratch/monoallelic/hjm31/bin_new/SNPsplit-0.5.0/SNPsplit \
--paired --bisulfite ./male_37 \
--snp_file ../snps/m37_snpsplit.txt \
trim_w37_1_bismark_bt2_pe.deduplicated.bam
 
 /scratch/monoallelic/hjm31/bin_new/SNPsplit-0.5.0/SNPsplit \
--paired --bisulfite -o ./queen_37 \
--snp_file ../snps/q37_snpsplit.txt \
trim_w37_1_bismark_bt2_pe.deduplicated.bam
 
 
 
 