import matplotlib.pyplot as plt
import numpy as np
import sys
import os

myfolder = sys.argv[1]

base_folder = '/home/lfrechet/BigBoy/capsid-assembly/HBV_KMC/%s/' % myfolder

sizes = []

nseed=100
for i in range(1,nseed+1):
    myfile = base_folder + 'seed=%d/prod/energy.dat' % i
    data = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols='NE')
    sizes.append(data['NE'][-1])
print(sizes)

max_size=250
mybins = np.arange(0-0.5,max_size+2-0.5,1)
hist, bins=np.histogram(sizes, bins=mybins)
bins = (bins[:-1]+bins[1:])/2
print('bins:',bins)

np.savetxt(base_folder + 'size_hist.txt', np.c_[bins,hist])

fig = plt.figure()
#plt.plot(bins, hist)
plt.bar(bins, hist,width=1)
#plt.xlim([75,130])
plt.xlim([100,500])
#plt.ylim([0,100])
plt.ylim([0,30])
plt.ylabel('counts')
plt.xlabel('assembly size (no. of dimers)')
plt.savefig('plots/size_hist_%s.png' % myfolder)
    
