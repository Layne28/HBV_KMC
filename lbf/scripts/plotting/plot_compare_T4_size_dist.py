import matplotlib.pyplot as plt
import numpy as np
import sys
import os

folders = ['T4_std', 'T4_std_lower_dg_conf']
#bins_list = [50,50,50,50]
bins_list = [30,30]
colors=['red','blue']
labels=[r'$\Delta g_{\text{conf}}=4.5$',r'$\Delta g_{\text{conf}}=4.0$' ]

fig = plt.figure()

for i in range(len(folders)):
    myfolder = folders[i]

    base_folder = '/home/lfrechet/BigBoy/capsid-assembly/HBV_KMC/%s/' % myfolder

    sizes = []

    nseed=100
    for j in range(1,nseed+1):
        myfile = base_folder + 'seed=%d/prod/energy.dat' % j
        data = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols='NE')
        sizes.append(data['NE'][-1])

    max_size=150
    mybins = np.arange(0-0.5,max_size+2-0.5,1)
    hist, bins=np.histogram(sizes, bins=mybins)
    #hist, bins=np.histogram(sizes, bins=bins_list[i])
    bins = (bins[:-1]+bins[1:])/2

    #np.savetxt(base_folder + 'size_hist.txt', np.c_[bins,hist])

    #plt.plot(bins, hist, color=colors[i], label=labels[i])
    plt.bar(bins, hist, color=colors[i], label=labels[i], alpha=0.5, width=1)
#plt.xlim([75,130])
plt.xlim([0,130])
#plt.ylim([0,100])
plt.ylim([0,100])
plt.legend()
plt.ylabel('counts')
plt.xlabel('assembly size (no. of dimers)')
plt.savefig('plots/T4_size_hists.png')
plt.show()
    
