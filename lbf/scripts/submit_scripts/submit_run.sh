#!/bin/bash

#SBATCH -J HBV_KMC          # Job name
#SBATCH -o HBV_KMC.o%j       # Name of stdout output file
#SBATCH -e HBV_KMC.e%j       # Name of stderr error file
#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 4
#SBATCH --mem-per-cpu=4G
#SBATCH -t 96:00:00        # Run time (hh:mm:ss)


#Run job
eval $1
