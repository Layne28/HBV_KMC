import numpy as np
import matplotlib.pyplot as plt
import sys
import os

base_folder = sys.argv[1]

print(os.listdir(base_folder))
subfolders = [e for e in os.listdir(base_folder) if e.startswith('seed')]
nseeds = len(subfolders)

fig = plt.figure()

for i in range(1,nseeds+1):
    #myfile = base_folder + '/seed=%d/energy.dat' % i
    myfile = base_folder + '/seed-%d/energy.dat' % i
    data = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols='NE') 
    print(data.shape)
    plt.plot(np.arange(data.shape[0]), data)
plt.savefig('avg_size_vs_time.png')
plt.show()