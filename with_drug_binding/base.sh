#!/bin/bash

#SBATCH --job-name=$JNAME
#SBATCH --output=$JNAME.txt

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --qos=medium

#SBATCH --time=20:00:00
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
./source/assemble $SEED 4200.000 40.000 800.000 0.240 0.480 $GB -11.5 0.02 -4.5 0.100 $dmud 0.0000 $drate 0.3 0.1 -0.1 0 -0.8 -.95
sleep 1
rm -rf source
echo "End $SEED @ `date`"
