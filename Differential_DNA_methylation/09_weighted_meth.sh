#!/bin/bash

#PBS -N weighted_meth
#PBS -l walltime=05:00:00
#PBS -l vmem=200gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=5

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed 
module load R/3.6.1 

R --save -q -f script.R