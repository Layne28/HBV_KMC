#!/bin/bash

#SBATCH -J HBV_KMC_dG       # Job name (do_submit_dG_sweep.sh overrides it per dG value)
#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH -N 1                # Total # of nodes
#SBATCH -n 1                # hbv_kmc is single-threaded
#SBATCH --mem-per-cpu=2G
#SBATCH -t 12:00:00         # Run time (hh:mm:ss); 1M sweeps take ~15 min on a workstation

# One seed of one dG_apoholo value. Submitted as a job array by
# do_submit_dG_sweep.sh, from the scratch directory that holds hbv_kmc,
# seeds.txt and the input file. The array index is the seed index: the line
# of seeds.txt to use, and the seed=<n> output subdirectory.

if [ "$#" -ne 1 ]; then
  echo "Usage: sbatch --array=1-<nseeds> $0 <input file>"
  exit 1
fi
infile=$1

./hbv_kmc ${infile} ${SLURM_ARRAY_TASK_ID} seeds.txt
