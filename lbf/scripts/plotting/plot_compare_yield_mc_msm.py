import pickle
import sys
import matplotlib.pyplot as plt
import os
import numpy as np

msm_file = sys.argv[1] #e.g. msm_tau=1.000000.pkl
tau = float((msm_file.split('.pkl')[0]).split('tau=')[-1])
print('tau: ', tau)
tau0=10000.0

#first get MC yield
times_to_capsid = []
myfolder = sys.argv[2]
subfolders = [e for e in os.listdir(myfolder) if e.startswith('seed')]
nseeds = len(subfolders)
print(myfolder)
for i in range(nseeds):
    myfile = myfolder + '/seed=%d/prod/energy.dat' % (i+1)
    if os.path.exists(myfile):
        #print(myfile)
        if os.path.getsize(myfile) > 0:
            print(myfile)
            data = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols=['sweep','NE'])

            #Get index where capsid first appears
            mask = (data['NE'] == 120)
            print(data['NE'])
            if mask.size!=0:
                first_occurrence_index = mask.argmax()
                if first_occurrence_index==0:
                    first_occurrence_index=-1
            else:
                first_occurrence_index = -1
            times_to_capsid.append(first_occurrence_index)


with open(msm_file, 'rb') as f:
    msm = pickle.load(f)

#print(msm.__dict__)
#print(msm.compute_committor([0],[20]))
MM = msm.__dict__['_MSM__macrostate_map']
#print(MM.__dict__)

MM_dict = MM.__dict__['_MacrostateMap__toIndex']

#print(MM_dict['ndimer=3_nCD=3'])

#print(MM.__dict__['ndimer=3_nCD=1'])

#Solve for time evolution
p0 = [0]*msm._MSM__num_states#[1,0] #initial distribution. 100% A
p0[0] = 1
T  = 20000  #final time (in lags) 
p  = msm.solve_FKE(p0, T)

#Get MC yield
mc_yield = np.zeros(T)
print(times_to_capsid)
for index in times_to_capsid:
    if index!=-1:
        print(index)
        dat_arr = np.zeros(T)
        dat_arr[index:] = 1
        mc_yield += dat_arr
mc_yield *= 1.0/100

#Plot time evolution of trimers and capsids
fig = plt.figure()
plt.plot(np.array(range(T+1))*(tau/tau0),p[:,MM_dict['ndimer=3_nCD=1']]+p[:,MM_dict['ndimer=3_nCD=3']],color='blue',label='trimer (MSM)')
plt.plot(np.array(range(T+1))*(tau/tau0),p[:,MM_dict['ndimer=120_nCD=60']],color='red',label='T=4 capsid (MSM)')
plt.plot(range(T),mc_yield,color='orange',linestyle='--',label='T=4 capsid (MC)')
plt.xlim([0,T])
plt.xlabel(r'MC sweeps/$10^4$')
plt.ylabel(r'probability')
plt.legend(loc='lower right')
plt.savefig('plots/T4_std_yield_vs_time_mc_vs_msm_tau=%f.png' % tau)
