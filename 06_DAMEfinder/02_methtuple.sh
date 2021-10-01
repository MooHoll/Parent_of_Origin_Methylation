#!/bin/bash

#PBS -N methtuple
#PBS -l walltime=16:00:00
#PBS -l vmem=40gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=2

# Run script in the working directory it was submitted in (56rs)
cd $PBS_O_WORKDIR 

# Make the methtuple environment
#module load python/gcc/2.7.13
#virtualenv --no-site-packages methtuple_env
#source /scratch/monoallelic/hjm32/bin/methtuple_env/bin/activate
# Download the github repository and copy to the methtuple_env
# Inside the new methtuple github folder run:
#python setup.py install

# Activate python enviroment
module load python/gcc/2.7.13
source /scratch/monoallelic/hjm32/bin/methtuple_env/bin/activate

for file in $(ls *methtuple.bam)
do
    methtuple --sc --gzip --mt CG -m 2 ${file}
done
