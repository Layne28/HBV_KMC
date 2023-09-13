#!/bin/bash

#SBATCH --job-name=$JNAME
#SBATCH --output=$JNAME.txt

#SBATCH --account=$ACCOUNT
#SBATCH --partition=$PARTITION
#SBATCH --qos=medium

#SBATCH --time=20:00:00
#SBATCH --exclude=compute-2-16,compute-0-25,compute-1-22,compute-2-19,compute-3-13,compute-2-14,compute-2-14,compute-2-16,compute-2-19,compute-2-10,compute-2-15,compute-2-18,compute-2-21,compute-2-23,compute-1-9,compute-1-18,compute-1-19,compute-3-13,compute-3-15,compute-3-16,compute-3-17,compute-3-21,compute-2-25


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
./source/assemble $SEED 4200.000 40.000 800.000 0.240 0.480 $GB $MU 0.02 $DMU 0.100 0.0 0.000 0.00000 0.3 0.1 -0.1 0 -0.8 -.95
sleep 1
rm -rf source
echo "End $SEED @ `date`"
