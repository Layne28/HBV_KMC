import matplotlib.pyplot as plt
import numpy as np
import sys
import os

folders = ['T4_std', 'CAM_high_salt','CAM_mod_salt','CAM_low_salt']
#bins_list = [50,50,50,50]
bins_list = [30,30,30,30]
colors=['black','green','red','blue']
labels=['no CAM, mod. salt','CAM, high salt', 'CAM, mod. salt','CAM, low salt']

fig = plt.figure()

for i in range(len(folders)):
    myfolder = folders[i]

    print(myfolder)

    base_folder = '/home/lfrechet/BigBoy/capsid-assembly/HBV_KMC/%s/' % myfolder

    sizes = []

    nseed=100
    for j in range(1,nseed+1):
        myfile = base_folder + 'seed=%d/prod/energy.dat' % j
        print(myfile)
        with open(myfile, 'r') as f:
            line = f.readline()
        print(line)
        data = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols='NE')
        sizes.append(data['NE'][-1])

    #max_size=2000
    #mybins = np.arange(0-0.5,max_size+2-0.5,1)
    #hist, bins=np.histogram(sizes, bins=mybins)
    hist, bins=np.histogram(sizes, bins=bins_list[i])
    bins = (bins[:-1]+bins[1:])/2

    #np.savetxt(base_folder + 'size_hist.txt', np.c_[bins,hist])

    plt.plot(bins, hist, color=colors[i], label=labels[i])
    #plt.bar(bins, hist,width=1)
#plt.xlim([75,130])
plt.xlim([0,600])
#plt.ylim([0,100])
plt.ylim([0,60])
plt.legend()
plt.ylabel('counts')
plt.xlabel('assembly size (no. of dimers)')
plt.savefig('plots/size_hists.png')
plt.show()
    
