'''
Get final size distribution
'''

import sys
import os
import numpy as np
import matplotlib.pyplot as plt 

#myfile = sys.argv[1] #energy.dat
base_folder = sys.argv[1] #should contain "seed" subfolders which in turn contain energy.dat

kadd_list = ['0.002','0.005','0.01','0.02']

plt.figure()

for kadd in kadd_list:

    print(float(kadd))
    frac = float(kadd_list[0])/float(kadd)

    sizes = []

    myfolder = base_folder + kadd + '/'

    subfolders = [e for e in os.listdir(myfolder) if e.startswith('seed')]
    nseeds = len(subfolders)

    for i in range(1,nseeds+1):
        myfile = myfolder + '/seed=%d/energy.dat' % (i)

        data = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols=('sweep','NE'))

        #final_size = data['NE'][-1]
        #mysize = data['NE'].shape[0]
        mysize = 360
        print(mysize)
        #print(frac)
        #print(int(mysize*frac))
        final_size = data['NE'][int(mysize*frac)-1]
        final_sweep = data['sweep'][-1]

        sizes.append(final_size)

    hist, bins = np.histogram(sizes, bins=np.arange(3,65))
    #hist, bins = np.histogram(sizes, bins=20)
    bin_centers = (bins[:-1] + bins[1:])/2
    print(bin_centers)

    mywidth=1.0
    plt.bar(bin_centers, hist,label='kadd=%s' % kadd,width=mywidth,alpha=0.5)

plt.legend()
plt.xlabel('size')
plt.ylabel('count')
plt.savefig('size_dist_lowsalt_vs_kadd.png')
plt.show()


        #print('seed: %d final sweep: %d and size: %d' % (i,final_sweep, final_size))

