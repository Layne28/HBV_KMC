#!/bin/bash

#SBATCH -J CG_KMC          # Job name
#SBATCH -o CG_KMC.o%j       # Name of stdout output file
#SBATCH -e CG_KMC.e%j       # Name of stderr error file
#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute-short             # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 1
#SBATCH -t 12:00:00        # Run time (hh:mm:ss)


if [ "$#" -eq 0 ]; then
  echo "Error: No arguments provided."
  echo "Usage: $0 <subfolder (e.g. T4_std, CAM_mod_salt)>"
  exit 1 # Exit with a non-zero status to indicate an error
fi
subfolder=$1

datetime=$(date +"%Y_%m_%d_%H_%M_%S")
scratch_dir="$SCRATCH/capsid-assembly/HBV_KMC/"
scratch_subdir="$SCRATCH/capsid-assembly/HBV_KMC/${subfolder}_analysis_$datetime/"
scratch_data_subdir="$SCRATCH/capsid-assembly/HBV_KMC/${subfolder}/"

#Make directories and copy files to SCRATCH
mkdir -p "${scratch_subdir}scripts/"
cp -r scripts/analysis ${scratch_subdir}scripts/
#cp input_files/${base_name}.in ${scratch_subdir}

#cd into SCRATCH directory
cd ${scratch_subdir}

#Launch jobs
python ${scratch_subdir}scripts/analysis/hbv_cg_traj_ABCD_sheet.py ${scratch_data_subdir}
