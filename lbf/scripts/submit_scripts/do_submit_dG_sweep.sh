#!/bin/bash

# Sweep dG_apoholo (G_holo - G_apo, in kT; used to model the oxidant
# concentration) for one base input file, with several seeds per value.
# Run from the lbf directory, like do_submit_run.sh:
#   scripts/submit_scripts/do_submit_dG_sweep.sh input_files/disassembly_crosslinked_extreme.in [nseeds]
# For each value it writes <base>_gox=<dG>.in in the scratch directory, with
# dG_apoholo set and output_dir pointing to $SCRATCH/capsid-assembly/HBV_KMC/<base>_gox=<dG>/,
# and submits one job array with one task per seed (see submit_dG_run.sh).

# dG_apoholo values (kT). In the equilibrated closed shell, switching a dimer
# to holo costs 44-80 kT for AB dimers and 107-161 kT for CD dimers, on top
# of dG_apoholo (stretch at the holo length plus dgintra). So:
#   -30            nothing switches (control)
#   -40 to -80     AB dimers start to switch, from the most stretched to most of them
#   -100 to -160   CD dimers start to switch too
dGs=(-30 -40 -45 -50 -55 -60 -70 -80 -100 -120 -140 -160)

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <base input file> [nseeds (default 10, at most the lines in seeds.txt)]"
  exit 1
fi
base_input=$1
nseeds=${2:-10}
if [ "$nseeds" -gt "$(wc -l < seeds.txt)" ]; then
  echo "Error: seeds.txt has only $(wc -l < seeds.txt) seeds"
  exit 1
fi

base_file=$(basename -- "$base_input")
base_name="${base_file%.*}"
base_name="${base_name%_gox=*}"   # an input already named for one value is treated as the base

scratch_dir="$SCRATCH/capsid-assembly/HBV_KMC/"
scratch_subdir="${scratch_dir}${base_name}_dG_sweep_scratch/"

# Copy the executable, seeds and scripts to SCRATCH
mkdir -p "${scratch_subdir}scripts/submit_scripts"
cp -r scripts/* ${scratch_subdir}scripts/
cp bin/hbv_kmc seeds.txt ${scratch_subdir}
cd ${scratch_subdir}

for dG in ${dGs[@]}; do
  name="${base_name}_gox=${dG}"
  infile="${name}.in"
  cp "${OLDPWD}/${base_input}" ${infile}
  sed -i "s|^output_dir:.*|output_dir: ${scratch_dir}${name}/|" ${infile}
  if grep -q "^dG_apoholo:" ${infile}; then
    sed -i "s|^dG_apoholo:.*|dG_apoholo: ${dG}|" ${infile}
  else
    echo "dG_apoholo: ${dG}" >> ${infile}
  fi
  echo "Submitting ${name}: seeds 1-${nseeds}"
  sbatch --array=1-${nseeds} --job-name="${name}" \
         --output="${name}_%a.out" --error="${name}_%a.err" \
         scripts/submit_scripts/submit_dG_run.sh ${infile}
done
