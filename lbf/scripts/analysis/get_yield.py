import pickle
import sys
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib as mpl

fig = plt.figure()
cmap = mpl.colormaps['plasma']
#colors = cmap(np.linspace(0,1,len(taus)+1))

#Get MC yield

times_to_T4_capsid = []
times_to_T3_capsid = []
myfolder = sys.argv[1]
subfolders = [e for e in os.listdir(myfolder) if e.startswith('seed')]
nseeds = len(subfolders)
print(myfolder)

skip_seeds = [41, 57]
for i in range(nseeds):
    myfile = myfolder + '/seed=%d/prod/energy.dat' % (i+1)
    if os.path.exists(myfile and (i+1) not in skip_seeds):
        #print(myfile)
        if os.path.getsize(myfile) > 0:
            print(myfile)
            data = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols=['sweep','NE','NCD_T4','NCD_T3'])

            #Get index where T4 capsid first appears
            mask_T4 = (data['NE'] == 120) & (data['NCD_T4']== 60)
            print(data['NE'])
            if mask_T4.size!=0:
                first_occurrence_index = mask_T4.argmax()
                if first_occurrence_index==0:
                    first_occurrence_index=-1
            else:
                first_occurrence_index = -1
            print(first_occurrence_index)
            times_to_T4_capsid.append(first_occurrence_index)

            #Get index where T3 capsid first appears
            mask_T3 = (data['NE'] == 90) & (data['NCD_T3']== 30)
            print(data['NE'])
            if mask_T3.size!=0:
                first_occurrence_index = mask_T3.argmax()
                if first_occurrence_index==0:
                    first_occurrence_index=-1
            else:
                first_occurrence_index = -1
            times_to_T3_capsid.append(first_occurrence_index)

#Get MC yield
tf = 2000000  #final time (/tau_0)
T4_yield = np.zeros(tf)
T3_yield = np.zeros(tf)
print(times_to_T4_capsid)
for index in times_to_T4_capsid:
    if index!=-1:
        print(index)
        dat_arr = np.zeros(tf)
        dat_arr[index:] = 1
        T4_yield += dat_arr
T4_yield *= 1.0/len(times_to_T4_capsid)

for index in times_to_T3_capsid:
    if index!=-1:
        print(index)
        dat_arr = np.zeros(tf)
        dat_arr[index:] = 1
        T3_yield += dat_arr
T3_yield *= 1.0/len(times_to_T3_capsid)

#Save to file
np.savetxt(myfolder + '/yield.txt', np.c_[np.array(range(tf)),T4_yield,T3_yield], header='time T4_yield T3_yield')

#Plot time evolution of capsids
plt.plot(range(tf),T4_yield,color='blue',linestyle='--',label='T4 yield')
plt.plot(range(tf),T3_yield,color='red',linestyle='--',label='T3 yield')
#plt.xlim([0,tf])
plt.ylim([0,1])
#plt.xscale('log')
#plt.xlim([100,3*10**4])
plt.xlabel(r'MC sweeps/$\tau_0$')
plt.ylabel(r'probability')
plt.legend(loc='lower right',fontsize=8,ncol=1)
plt.savefig('plots/MC_T4_T3_yield.png')
plt.show()
