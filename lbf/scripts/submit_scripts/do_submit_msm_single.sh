#!/bin/bash

folders=('mod_salt_lower_gdc_high_freq')
lag_times=(10 20 50 100 200) #(200 800 1400 2000)

do_pad=1
do_abs=0

#Submit jobs
for folder in ${folders[@]}; do
    for lag_time in ${lag_times[@]}; do
        echo "Submitting job: ${folder} ${lag_time}"
        sbatch scripts/submit_scripts/submit_msm.sh ${folder} ${lag_time} ${do_pad} ${do_abs}
    done
done

