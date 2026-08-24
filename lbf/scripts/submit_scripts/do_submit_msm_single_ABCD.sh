#!/bin/bash

folders=('CAM_mod_salt')
lag_times=(400 600 1000 1200 1600 1800) #(200 800 1400 2000) #(2000) #(200 800 1400 2000)

do_pad=1
do_abs=0
target_state='none'

#Submit jobs
for folder in ${folders[@]}; do
    for lag_time in ${lag_times[@]}; do
        echo "Submitting job: ${folder} ${lag_time}"
        sbatch scripts/submit_scripts/submit_msm_ABCD.sh ${folder} ${lag_time} ${do_pad} ${do_abs} ${target_state}
    done
done

