# Issue: the alignments to the alternate reference genomes didn't pull through the chromosome name,
# just a corresponding number, this means the output of the differential meth have no chromosome information
# its lazy, but as there are only a handful I will just look them up manually 

# Get chromosome names and corresponding numbers from an alternate reference genome 
grep ">" m08_unique.fasta > chrom_numbers.txt

sed -i 's/>//g' chrom_numbers.txt

# Pull out each position and find the corresponding chromosome name etc.
grep "14" chrom_numbers.txt | head


# 14 = NC_015775.1
# 179 = NW_003565564.1
# 18 = NC_015779.1
# 50  = NW_003565435.1
# 7 = NC_015768.1
# 3 = NC_015764.1
# 5 = NC_015766.1
# 58 = NW_003565443.1

# Replace the chromosome number with the name in the methylkit output file
sed -i 's/,14,/,NC_015775.1,/' DMRs_min5percentDiff_qval0.05_MSCfilter.csv
sed -i 's/,179,/,NW_003565564.1,/g' DMRs_min5percentDiff_qval0.05_MSCfilter.csv
sed -i 's/,18,/,NC_015779.1,/' DMRs_min5percentDiff_qval0.05_MSCfilter.csv
sed -i 's/,50,/,NW_003565435.1,/' DMRs_min5percentDiff_qval0.05_MSCfilter.csv
sed -i 's/,7,/,NC_015768.1,/' DMRs_min5percentDiff_qval0.05_MSCfilter.csv
sed -i 's/,3,/,NC_015764.1,/' DMRs_min5percentDiff_qval0.05_MSCfilter.csv
sed -i 's/,5,/,NC_015766.1,/' DMRs_min5percentDiff_qval0.05_MSCfilter.csv
sed -i 's/,58,/,NW_003565443.1,/' DMRs_min5percentDiff_qval0.05_MSCfilter.csv


### LATER:
# Needed to change chromosome names for whole file so...

sed -i 's/:1//g' chrom_numbers.txt