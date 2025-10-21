import matplotlib.pyplot as plt
import numpy as np
import sys

second_col='NCD_T4'

#myfile = base_folder + '/seed-%d/energy.dat' % i
myfile = '/home/lfrechet/BigBoy/capsid-assembly/HBV_KMC/mod_salt/seed=9/prod/energy.dat'
data = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols=('NE',second_col))
data_CD =  np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols=('NCD_Hex','NCD_other','NCD_T4','NCD_T3'))
numerical_data = np.array([data['NE'], data[second_col]]).T
#numerical_data_CD = np.array([data_CD['NCD_Hex'], data_CD['NCD_other'], data_CD['NCD_T4'], data_CD['NCD_T3']]).T
numerical_data_CD = np.array([data_CD['NCD_Hex'], data_CD['NCD_other'], data_CD['NCD_T4']]).T
numerical_data_CD_tot = np.sum(numerical_data_CD, axis=1)
print(numerical_data_CD)
print(numerical_data.shape)
print(numerical_data_CD.shape)

fig = plt.figure()
plt.scatter(numerical_data[:,0],numerical_data_CD_tot[:])
plt.scatter(numerical_data[:,0],numerical_data[:,1])
plt.show()