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

max_sweep =0

for i in range(1,nseeds+1):
    myfile = base_folder + 'seed=%d/energy.dat' % i

    data = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols=('sweep','NE'))

    final_size = data['NE'][-1]
    final_sweep = data['sweep'][-1]

    if final_sweep > max_sweep:
        max_sweep = final_sweep

    print('seed: %d final sweep: %d and size: %d' % (i,final_sweep, final_size))

data_list = []
sweep_data = []
for i in range(1, nseeds+1):
    myfile = base_folder + 'seed=%d/energy.dat' % i

    data = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols=('sweep','NE'))

    final_size = data['NE'][-1]
    final_sweep = data['sweep'][-1]

    if final_sweep==max_sweep:
        sweep_data = data['sweep']
        #Remove duplicate rows
        select_data = np.c_[data['sweep'],data['NE']]
        select_data = select_data[np.unique(select_data[:,0],return_index=True,axis=0)[1]]
        #select_data = np.unique(select_data,axis=0)
        print(i, select_data.shape)
        if select_data.shape[0]==201:
            data_list.append(select_data[:,1])

print(len(data_list))
avg_size = sum(data_list)/len(data_list)

np.savetxt(base_folder + '/size_vs_time.txt', np.c_[sweep_data,avg_size])

fig = plt.figure()
plt.plot(sweep_data, avg_size,linewidth=2,color='black')
for data in data_list:
    plt.plot(sweep_data, data,linewidth=0.5)
plt.savefig('test.png')
#plt.show()


