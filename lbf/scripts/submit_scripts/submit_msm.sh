#!/bin/bash

#SBATCH -J MSM_KMC          # Job name
#SBATCH -o MSM_KMC.o%j       # Name of stdout output file
#SBATCH -e MSM_KMC.e%j       # Name of stderr error file
#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute-short             # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 5
#SBATCH --mem-per-cpu=2G
#SBATCH -t 12:00:00        # Run time (hh:mm:ss)


if [ "$#" -eq 0 ]; then
  echo "Error: No arguments provided."
  echo "Usage: $0 <subfolder (e.g. T4_std, CAM_mod_salt)> <ndimer interval> <nCD interval> <ndrug interval> <lag time>"
  exit 1 # Exit with a non-zero status to indicate an error
fi
subfolder=$1
lag=$2
do_pad=$3
do_abs=$4

datetime=$(date +"%Y_%m_%d_%H_%M_%S")
scratch_dir="$SCRATCH/capsid-assembly/HBV_KMC/"
scratch_subdir="$SCRATCH/capsid-assembly/HBV_KMC/${subfolder}_tau=${lag}_analysis_$datetime/"
scratch_data_subdir="$SCRATCH/capsid-assembly/HBV_KMC/${subfolder}/"

#Make directories and copy files to SCRATCH
mkdir -p "${scratch_subdir}scripts/"
cp -r scripts/analysis ${scratch_subdir}scripts/
#cp input_files/${base_name}.in ${scratch_subdir}

#cd into SCRATCH directory
cd ${scratch_subdir}

#Launch jobs
python ${scratch_subdir}scripts/analysis/hbv_msm.py ${scratch_data_subdir} ${lag} ${do_pad} ${do_abs}
