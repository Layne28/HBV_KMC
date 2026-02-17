#!/bin/bash

#command line arguent: name of base input file (e.g. T4.in)
input_base_file=$1 

filename=$(basename -- "$input_base_file")
extension="${filename##*.}"
base_name="${filename%.*}"

echo $filename

startseed=$2
endseed=$3
seeds=($(seq $startseed $endseed))
#mudrugs=(-1.0)

mylen=${#seeds[@]}
echo $mylen

cmd_list=()

for seed in "${seeds[@]}"
do
    mycommand="./hbv_kmc ${filename} $seed many_seeds.txt"
    cmd_list+=("${mycommand}")
done

rm -f "${base_name}_ManySweeps.txt"
for cmd in "${cmd_list[@]}"
do
    echo "${cmd}"
    echo ""
    printf "${cmd}\n" >> "${base_name}_ManySweeps.txt"
done
