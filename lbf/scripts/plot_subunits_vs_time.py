import numpy as np
import matplotlib.pyplot as plt

base_folder = '/home/lfrechet/research/capsid-assembly/hbv_assembly/disassembly/HBV_KMC/lbf/data/'
subfolders = ['disassembly_crosslinked_extreme_kl=0.02','disassembly_crosslinked_extreme','disassembly_crosslinked_extreme_kl=0.0002']
colors=['red','purple','blue']
labels=[r'$k_l^0=2\times 10^{-2}$',r'$k_l^0=2\times 10^{-3}$',r'$k_l^0=2\times 10^{-4}$']

fig = plt.figure()
for s, c, l in zip(subfolders, colors, labels):
    data = np.loadtxt(base_folder + s + '/seed=1/prod/energy.dat', delimiter=',', skiprows=1)
    plt.plot(data[:,0],data[:,17],color=c,label=l)
plt.xlabel('MC sweeps')
plt.ylabel('No. of dimers')
plt.legend(loc='lower left')
plt.ylim([50,121])
plt.xlim([-3000,400000])
plt.savefig('plots/disassembly_no_dimers_vs_time.png')