import pickle
import sys
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib as mpl
sys.path.insert(1, '/home/lfrechet/research/capsid-assembly/hbv_assembly/HBV_KMC/lbf/scripts/analysis')
import MSM

msm_file = sys.argv[1] #msm_tau=XXX.pkl
with open(msm_file, 'rb') as f:
        msm = pickle.load(f)

tau0 = 100.0
tau = float((msm_file.split('.pkl')[0]).split('tau=')[-1])
print('tau: ', tau)
tf = 2000000#2000000  #final time (/tau0)

MM = msm.__dict__['_MSM__macrostate_map']
MM_dict = MM.__dict__['_MacrostateMap__toIndex']

#Solve for time evolution
p0 = [0]*msm._MSM__num_states#[1,0] #initial distribution. 100% A
print('loc: ', MM_dict['ndimer=3_nCD=1'])
p0[MM_dict['ndimer=3_nCD=1']] = 1 #check this
T = int(tf/(tau/tau0))
p  = msm.solve_FKE(p0, T)

print(p[-1,:])
print(p[-1,:].argmax(), p[-1,:].max())
print(MM_dict['ndimer=120_nCD=60'])
print(MM_dict)

fig = plt.figure()
cmap = mpl.colormaps['plasma']
#colors = cmap(np.linspace(0,1,len(taus)+1))

#Plot time evolution of trimers and capsids
plt.plot(np.array(range(T+1))*(tau/tau0),p[:,MM_dict['ndimer=120_nCD=60']],color='red',label=r'MSM, $\tau/\tau_0=%.0f$' % (tau/tau0))
#plt.xlim([0,tf])
#plt.ylim([0,1])
#plt.xscale('log')
#plt.xlim([100,3*10**4])
plt.xlabel(r'MC sweeps/$\tau_0$')
plt.ylabel(r'probability')
plt.legend(loc='lower right',fontsize=8,ncol=1)
plt.savefig('plots/msm_states_vs_time.png')
#plt.show()
