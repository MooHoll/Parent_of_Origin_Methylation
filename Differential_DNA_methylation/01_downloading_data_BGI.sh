############################################################################
### Getting data from BGI ftp servers
############################################################################

wget <ftp download link for file or directory> --ftp-user=<username> --ftp-password=<password>

ftp://cdts-wh.genomics.cn/Clean_Data/10895-1/* 
ftp://cdts-wh.genomics.cn/Clean_Data/10895-10/*
ftp://cdts-wh.genomics.cn/Clean_Data/10895-11/*
ftp://cdts-wh.genomics.cn/Clean_Data/10895-12/*
ftp://cdts-wh.genomics.cn/Clean_Data/10895-2/*
ftp://cdts-wh.genomics.cn/Clean_Data/10895-3/*
ftp://cdts-wh.genomics.cn/Clean_Data/10895-4/*
ftp://cdts-wh.genomics.cn/Clean_Data/10895-5/*
ftp://cdts-wh.genomics.cn/Clean_Data/10895-6/*
ftp://cdts-wh.genomics.cn/Clean_Data/10895-7/*
ftp://cdts-wh.genomics.cn/Clean_Data/10895-8/*
ftp://cdts-wh.genomics.cn/Clean_Data/10895-9/*

# Re-sequenced ones:
ftp://cdts-wh.genomics.cn/Clean_Data/new/10895-1/*
ftp://cdts-wh.genomics.cn/Clean_Data/new/10895-8/*


############################################################################

#!/bin/bash

#PBS -N downloading_data
#PBS -l walltime=60:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

wget -r ftp://cdts-wh.genomics.cn/Clean_Data/10895-1/* --ftp-user=20180327F18FTSEUHT0130 --ftp-password=BEEmjdM
wget -r ftp://cdts-wh.genomics.cn/Clean_Data/10895-10/* --ftp-user=20180327F18FTSEUHT0130 --ftp-password=BEEmjdM
wget -r ftp://cdts-wh.genomics.cn/Clean_Data/10895-11/* --ftp-user=20180327F18FTSEUHT0130 --ftp-password=BEEmjdM
wget -r ftp://cdts-wh.genomics.cn/Clean_Data/10895-12/* --ftp-user=20180327F18FTSEUHT0130 --ftp-password=BEEmjdM
wget -r ftp://cdts-wh.genomics.cn/Clean_Data/10895-2/* --ftp-user=20180327F18FTSEUHT0130 --ftp-password=BEEmjdM
wget -r ftp://cdts-wh.genomics.cn/Clean_Data/10895-3/* --ftp-user=20180327F18FTSEUHT0130 --ftp-password=BEEmjdM
wget -r ftp://cdts-wh.genomics.cn/Clean_Data/10895-4/* --ftp-user=20180327F18FTSEUHT0130 --ftp-password=BEEmjdM
wget -r ftp://cdts-wh.genomics.cn/Clean_Data/10895-5/* --ftp-user=20180327F18FTSEUHT0130 --ftp-password=BEEmjdM
wget -r ftp://cdts-wh.genomics.cn/Clean_Data/10895-6/* --ftp-user=20180327F18FTSEUHT0130 --ftp-password=BEEmjdM
wget -r ftp://cdts-wh.genomics.cn/Clean_Data/10895-7/* --ftp-user=20180327F18FTSEUHT0130 --ftp-password=BEEmjdM
wget -r ftp://cdts-wh.genomics.cn/Clean_Data/10895-8/* --ftp-user=20180327F18FTSEUHT0130 --ftp-password=BEEmjdM
wget -r ftp://cdts-wh.genomics.cn/Clean_Data/10895-9/* --ftp-user=20180327F18FTSEUHT0130 --ftp-password=BEEmjdM


############################################################################

#!/bin/bash

#PBS -N downloading_data
#PBS -l walltime=00:30:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

wget -r ftp://cdts-wh.genomics.cn/Clean_Data/new/10895-1/* --ftp-user=20180327F18FTSEUHT0130 --ftp-password=BEEmjdM
wget -r ftp://cdts-wh.genomics.cn/Clean_Data/new/10895-8/* --ftp-user=20180327F18FTSEUHT0130 --ftp-password=BEEmjdM


############################################################################

# After all data is downloaded, need to gunzip any that were run on multiple lanes, use cat to make one file
# and then gzip again
