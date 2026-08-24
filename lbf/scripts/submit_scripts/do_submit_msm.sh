#!/bin/bash

folders=('high_salt' 'mod_salt' 'low_salt' 'CAM_high_salt' 'CAM_mod_salt' 'CAM_low_salt')
#lag_times=(200 400 600 800 1000 1200 1400 1600 1800 2000)
lag_times=(200 800 1400 2000)

do_pad=1
do_abs=1

#Submit jobs
for folder in ${folders[@]}; do
    for lag_time in ${lag_times[@]}; do
        echo "Submitting job: ${folder} ${lag_time}"
        sbatch scripts/submit_scripts/submit_msm.sh ${folder} ${lag_time} ${do_pad} ${do_abs}
    done
done

