#!/bin/bash

#command line arguent: (1) name of base input file (e.g. T4.in)
#####################  (2) name of scratch folder containing out files
input_base_file=$1 
input_scratch_folder=$2

filename=$(basename -- "$input_base_file")
extension="${filename##*.}"
base_name="${filename%.*}"

nseed=100
bad_string="Segmentation"
bad_string2="Aborted"
bad_string3="Killed"
bad_string4="Error"

echo $filename

cmd_list=()

for file in "${input_scratch_folder}HBV_KMC_"*"err"*; do
  #echo "Processing file: $file"

  # Get the last line of the file
  last_line=$(tail -n 1 "$file")
  if [[ "$last_line" =~ "$bad_string" || "$last_line" == *"$bad_string2"* || "$last_line" == *"$bad_string3"* || "$last_line" == *"$bad_string4"* ]]; then
    echo "Substring '$bad_string' or '$bad_string2' or '$bad_string3' or '$bad_string4' found in the last line of '$file'."
    #number=$(echo "$file" | sed -E 's/.*out([[:digit:]]+)/\1/')
    #myline=$(sed '4q;d' $file)
    number=$(echo "$file" | sed -E 's/.*HBV_KMC_([[:digit:]]+).err/\1/')
    #number=$(echo  #"${tmp%.err*}"
    #number=$((number + 1))
    echo "seed: $number"
    mycommand="./hbv_kmc ${filename} $number seeds.txt"
    cmd_list+=("${mycommand}")
    echo "Moving old data to subfolder."
    myfolder="/scratch0/laynefrechette/capsid-assembly/HBV_KMC/${base_name}/seed=${number}/"
    mkdir -p "${myfolder}/segfault_data"
    mv ${myfolder}/* ${myfolder}/segfault_data/

    mkdir -p ${input_scratch_folder}/segfault_data/
    mv ${file} ${input_scratch_folder}/segfault_data/
    #echo "moving ${myfolder}"
  fi
done

#Also go through and check that seed directories exist
for ((i=1; i < $((nseed+1)); i++)); do
    myfolder="/scratch0/laynefrechette/capsid-assembly/HBV_KMC/${base_name}/seed=$i/"
    if [ ! -d "${myfolder}/prod" ]; then
        echo "Directory '$myfolder' does not exist."
        mycommand="./hbv_kmc ${filename} $i seeds.txt"

        # Flag to track if the string is found
	found=0

	# Loop through the array to check if the string exists
	for element in "${cmd_list[@]}"; do
	  if [[ "$element" == "$mycommand" ]]; then
	    found=1
	    break # Exit the loop if found
	  fi
	done

	# If the string was not found, append it to the array
	if [[ "$found" -eq 0 ]]; then
	  cmd_list+=("$mycommand")
	  echo "Added '$mycommand' to the array."
	else
	  echo "'$mycommand' is already in the array."
	fi
    fi
done

rm -f "${base_name}_redoSweeps.txt"
for cmd in "${cmd_list[@]}"
do
    echo "${cmd}"
    echo ""
    printf "${cmd}\n" >> "${base_name}_redoSweeps.txt"
done
