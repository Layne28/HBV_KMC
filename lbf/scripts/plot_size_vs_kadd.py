import numpy as np
import matplotlib.pyplot as plt

kadd_list = ['0.002','0.005','0.01','0.02']

fig = plt.figure()
for kadd in kadd_list:
    data = np.loadtxt('/scratch/09029/tg883285/capsid-assembly/HBV_KMC/smriti_drug_low_salt_kadd=%s/size_vs_time.txt' % kadd)
    print((float(kadd)/float(kadd_list[-1])))
    plt.scatter(data[:,0]*(float(kadd)/float(kadd_list[-1])),data[:,1],label='kadd=%s' % kadd)
    #plt.plot(data[:,0],data[:,1],label='kadd=%s' % kadd)
plt.xlabel('MC sweeps')
plt.ylabel('average size')
plt.legend()
plt.savefig('size_vs_time_compare_kadd_scaled.png')
#plt.savefig('size_vs_time_compare_kadd.png')

