#!/bin/bash

folders=('mod_salt_farri_gdc')
lag_times=(200) #(200 800 1400 2000)

do_pad=0
do_abs=0

#Submit jobs
for folder in ${folders[@]}; do
    for lag_time in ${lag_times[@]}; do
        echo "Submitting job: ${folder} ${lag_time}"
        sbatch scripts/submit_scripts/submit_msm.sh ${folder} ${lag_time} ${do_pad} ${do_abs}
    done
done

