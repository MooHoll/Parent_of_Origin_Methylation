# On NCBI on list of datasets, send to, file, accession list, put list in nano on ALICE
#module load sratoolkit/2.9.6
#prefetch --option-file accessions.txt
# Lol of course this doesn't work, after a quick google it looks like the sratoolkit needs updating on ALICE

#  888 - Africanized mother x European father from block A. 
#  875 - European mother x Africanized father from block A. 
#  894 - Africanized mother x European father from block B. 
#  882 - European mother x Africanized father from block B.

# For now, get WGBS data (this is the raw data from the Xu et al. study)
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11148038/SRR11148038 #882_5
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11148039/SRR11148039 #882_3
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11148040/SRR11148040 #882_2
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11148041/SRR11148041 #882_0
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11148042/SRR11148042 #875_2
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11148043/SRR11148043 #875_1
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11148044/SRR11148044 #888_F
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11148045/SRR11148045 #888_E
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11148046/SRR11148046 #894_6
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11148047/SRR11148047 #894_5
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11148048/SRR11148048 #894_2
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11148049/SRR11148049 #894_1
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11148050/SRR11148050 #888_2
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11148051/SRR11148051 #888_1
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11148052/SRR11148052 #875_F
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11148053/SRR11148053 #875_E

# Get the genome data of the parents from the Galbraith et al. study
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR3037350/SRR3037350 #875 queen
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR3037351/SRR3037351 #875 drone
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR3037352/SRR3037352 #888 queen
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR3037353/SRR3037353 #888 drone
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR3037354/SRR3037354 #882 queen
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR3037355/SRR3037355 #882 drone
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR3037356/SRR3037356 #894 queen
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR3037357/SRR3037357 #894 drone

# Load software to convert SRR to fastq files
module load sratoolkit/2.9.6

# Convert to fastq file using the list of accessions we got earlier
cat sra_list.txt | while read line
do
    fasterq-dump --split-files ${line}
done

gzip SRR*