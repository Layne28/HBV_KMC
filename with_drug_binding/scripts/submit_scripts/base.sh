#!/bin/bash

#SBATCH --job-name=$JNAME
#SBATCH --output=$JNAME.txt

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --nodelist=compute-9-2
#SBATCH --qos=medium

#SBATCH --time=2-12:00:00
#SBATCH --ntasks=1
#SBATCH --chdir=$DIRF

#module load GSL/2.5
module load share_modules/GSL/2.5 
echo "Start $SEED @ `date`"
scp -r ../source $DIRX ;
sleep 10
cd $DIRX ;
cd source ;
make ;
sleep 10 
cd .. ;
# 0.404 0.247 -0.146 0 -0.663 -.8
./source/assemble $SEED 4200.000 40.000 800.000 0.240 0.480 $GB -12.2 0.02 $DMU 0.100 0.0 0.04 $kd 0.3 0.1 -0.1 0 -0.4 -.95
sleep 1
rm -rf source
echo "End $SEED @ `date`"
