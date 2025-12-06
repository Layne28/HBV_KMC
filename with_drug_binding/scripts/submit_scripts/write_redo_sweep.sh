#!/bin/bash

#command line arguent: (1) name of base input file (e.g. T4.in)
#####################  (2) name of scratch folder containing out files
input_base_file=$1 
input_scratch_folder=$2

filename=$(basename -- "$input_base_file")
extension="${filename##*.}"
base_name="${filename%.*}"

bad_string="Segmentation"

echo $filename

cmd_list=()

for file in "${input_scratch_folder}out"*; do
  #echo "Processing file: $file"

  # Get the last line of the file
  last_line=$(tail -n 1 "$file")
  if [[ "$last_line" =~ "$bad_string" ]]; then
    echo "Substring '$bad_string' found in the last line of '$file'."
    #number=$(echo "$file" | sed -E 's/.*out([[:digit:]]+)/\1/')
    myline=$(sed '4q;d' $file)
    number=$(echo "$myline" | sed -E 's/.*seed=([[:digit:]]+)/\1/')
    #number=$((number + 1))
    echo "seed: $number"
    mycommand="./hbv_kmc ${filename} $number seeds.txt"
    cmd_list+=("${mycommand}")
    echo "Moving old data to subfolder."
    myfolder="/scratch/09029/tg883285/capsid-assembly/HBV_KMC/${base_name}/seed=${number}/"
    mkdir -p "${myfolder}/segfault_data"
    mv ${myfolder}/* ${myfolder}/segfault_data/
    #echo "moving ${myfolder}"
  fi
done


rm -f "${base_name}_redoSweeps.txt"
for cmd in "${cmd_list[@]}"
do
    echo "${cmd}"
    echo ""
    printf "${cmd}\n" >> "${base_name}_redoSweeps.txt"
done
