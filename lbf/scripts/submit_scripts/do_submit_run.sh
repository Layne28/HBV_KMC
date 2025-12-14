#!/bin/bash

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
scratch_subdir="$SCRATCH/capsid-assembly/HBV_KMC/${base_name}_scratch/"
#$datetime/"

#Make directories and copy files to SCRATCH
mkdir -p "${scratch_subdir}scripts/submit_scripts"
cp -r scripts/* ${scratch_subdir}scripts/
cp ${sweeps_file} ${scratch_subdir}/
cp bin/* ${scratch_subdir}
cp seeds.txt ${scratch_subdir}
infile=$(head -n 1 ${sweeps_file} | awk '{print $2}')
cp input_files/$infile ${scratch_subdir}
#cp input_files/${base_name}.in ${scratch_subdir}

#cd into SCRATCH directory
cd ${scratch_subdir}

#modify .in file to write to scratch
base_name_no_redo="${base_name%_redo}"
new_out_folder="${scratch_dir}${base_name_no_redo}/"
sed -i "s|output_dir:.*|output_dir: ${new_out_folder}|" ${scratch_subdir}${base_name_no_redo}.in

#Submit jobs
while read l; do
    echo "Submitting job: $l"
    #Get seed
    seed=$(echo "$l" | cut -d " " -f 3)
    echo "seed: ${seed}"
    sbatch --output="HBV_KMC_${seed}.out" --error="HBV_KMC_${seed}.err" scripts/submit_scripts/submit_run.sh "$l"
done <"${sweeps_file}"

