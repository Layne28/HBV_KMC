#!/bin/bash

#SBATCH -J HBV_KMC          # Job name
#SBATCH -o HBV_KMC.o%j       # Name of stdout output file
#SBATCH -e HBV_KMC.e%j       # Name of stderr error file
#SBATCH -p spr,icx,skx             # Queue (partition) name
#SBATCH -N 4               # Total # of nodes 
##SBATCH -N 1               # Total # of nodes 
##SBATCH -t 24:00:00        # Run time (hh:mm:ss)
#SBATCH -t 12:00:00        # Run time (hh:mm:ss)
##SBATCH --mail-user=username@tacc.utexas.edu
##SBATCH --mail-type=all    # Send email at begin and end of job

module load pylauncher
LD_LIBRARY_PATH=$HOME/.local/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH


if [ "$#" -eq 0 ]; then
  echo "Error: No arguments provided."
  echo "Usage: $0 <sweeps file>"
  exit 1 # Exit with a non-zero status to indicate an error
fi
sweeps_file=$1
base_file=$(basename -- "$sweeps_file")
base_name="${base_file%.*}"
base_name="${base_name%Sweeps}"
echo $base_name

datetime=$(date +"%Y_%m_%d_%H_%M_%S")
scratch_dir="$SCRATCH/capsid-assembly/HBV_KMC/"
scratch_subdir="$SCRATCH/capsid-assembly/HBV_KMC/${base_name}_$datetime/"

#Make directories and copy files to SCRATCH
mkdir -p "${scratch_subdir}scripts/submit_scripts"
cp -r scripts/* ${scratch_subdir}scripts/
cp ${sweeps_file} ${scratch_subdir}/
cp bin/* ${scratch_subdir}
cp seeds.txt ${scratch_subdir}
cp -r "source" ${scratch_subdir}
infile=$(head -n 1 ${sweeps_file} | awk '{print $2}')
cp input_files/$infile ${scratch_subdir}
#cp input_files/${base_name}.in ${scratch_subdir}

#cd into SCRATCH directory
cd ${scratch_subdir}

#modify .in file to write to scratch
new_out_folder="${scratch_dir}${base_name}/"
sed -i "s|output_dir:.*|output_dir: ${new_out_folder}|" ${scratch_subdir}${base_name}.in

#Launch jobs
python ${scratch_subdir}scripts/submit_scripts/launch_sweep.py ${scratch_subdir}${sweeps_file}
