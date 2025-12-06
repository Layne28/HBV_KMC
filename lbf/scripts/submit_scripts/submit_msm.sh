#!/bin/bash

#SBATCH -J HBV_KMC          # Job name
#SBATCH -o HBV_KMC.o%j       # Name of stdout output file
#SBATCH -e HBV_KMC.e%j       # Name of stderr error file
#SBATCH -p spr,icx,skx             # Queue (partition) name
##SBATCH -N 4               # Total # of nodes 
#SBATCH -N 1               # Total # of nodes 
##SBATCH -t 24:00:00        # Run time (hh:mm:ss)
#SBATCH -t 1:00:00        # Run time (hh:mm:ss)
##SBATCH --mail-user=username@tacc.utexas.edu
##SBATCH --mail-type=all    # Send email at begin and end of job

module load pylauncher
LD_LIBRARY_PATH=$HOME/.local/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH


if [ "$#" -eq 0 ]; then
  echo "Error: No arguments provided."
  echo "Usage: $0 <subfolder (e.g. T4_std, CAM_mod_salt)> <ndimer interval> <nCD interval> <ndrug interval> <lag time>"
  exit 1 # Exit with a non-zero status to indicate an error
fi
subfolder=$1
ndimer=$2
nCD=$3
ndrug=$4
lag=$5

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
python ${scratch_subdir}scripts/analysis/hbv_msm.py ${scratch_data_subdir} ${ndimer} ${nCD} ${ndrug} ${lag}
