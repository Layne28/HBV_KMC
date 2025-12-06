'''
Get final size distribution
'''

import sys
import os
import numpy as np
import matplotlib.pyplot as plt 

#myfile = sys.argv[1] #energy.dat
base_folder = sys.argv[1] #should contain "seed" subfolders which in turn contain energy.dat

subfolders = [e for e in os.listdir(base_folder) if e.startswith('seed')]
nseeds = len(subfolders)

for i in range(1,nseeds+1):
    myfile = base_folder + 'seed=%d/energy.dat' % i

    data = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols=('sweep','NE'))

    final_size = data['NE'][-1]
    final_sweep = data['sweep'][-1]

    print('seed: %d final sweep: %d and size: %d' % (i,final_sweep, final_size))

