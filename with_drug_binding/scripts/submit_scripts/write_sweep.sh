#!/bin/bash

nseed=$1
seeds=($(seq 1 $nseed))

base_name="farri_T4"

mylen=${#seeds[@]}
echo $mylen

cmd_list=()

mapfile -t seedvalues < seeds.txt

counter=1
for seed in "${seedvalues[@]}"
do
    echo $seed
    mycommand="./source/assemble $seed 4200.000 40.000 800.000 0.0 0.480 -9.800 -11.500 0.0200 -4.500 0.100 1.000 1.000 0.000000 0.300 0.100 -0.100 0.000 -0.800 -0.950 seed=${counter}"
    cmd_list+=("${mycommand}")
    ((counter++))
done

rm -f "${base_name}Sweeps.txt"
for cmd in "${cmd_list[@]}"
do
    echo "${cmd}"
    echo ""
    printf "${cmd}\n" >> "${base_name}Sweeps.txt"
done
