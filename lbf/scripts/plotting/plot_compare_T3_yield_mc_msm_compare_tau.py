import pickle
import sys
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib as mpl

taus = [10000.0, 50000.0, 80000.0, 100000.0, 300000.0]

tf = 20000  #final time (/10^4)

fig = plt.figure()
cmap = mpl.colormaps['plasma']
colors = cmap(np.linspace(0,1,len(taus)+1))

#first get MC yield
times_to_capsid = []
myfolder = sys.argv[1]
subfolders = [e for e in os.listdir(myfolder) if e.startswith('seed')]
nseeds = len(subfolders)
print(myfolder)
for i in range(nseeds):
    myfile = myfolder + '/seed=%d/prod/energy.dat' % (i+1)
    if os.path.exists(myfile):
        #print(myfile)
        if os.path.getsize(myfile) > 0:
            print(myfile)
            data = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols=['sweep','NE','NCD_T3'])

            #Get index where capsid first appears
            mask = (data['NE'] == 90) & (data['NCD_T3']== 30)
            print(data['NE'])
            if mask.size!=0:
                first_occurrence_index = mask.argmax()
                if first_occurrence_index==0:
                    first_occurrence_index=-1
            else:
                first_occurrence_index = -1
            times_to_capsid.append(first_occurrence_index)

#Get MC yield
mc_yield = np.zeros(tf)
print(times_to_capsid)
for index in times_to_capsid:
    if index!=-1:
        print(index)
        dat_arr = np.zeros(tf)
        dat_arr[index:] = 1
        mc_yield += dat_arr
mc_yield *= 1.0/100

#Now collect MSM yields
tau0=10000.0
cnt = 0
for tau in taus:

    msm_file = myfolder + '/msm_tau=%f.pkl' % tau
    print('tau: ', tau)


    with open(msm_file, 'rb') as f:
        msm = pickle.load(f)

    #print(msm.__dict__)
    #print(msm.compute_committor([0],[20]))
    MM = msm.__dict__['_MSM__macrostate_map']
    MM_dict = MM.__dict__['_MacrostateMap__toIndex']

    #Solve for time evolution
    p0 = [0]*msm._MSM__num_states#[1,0] #initial distribution. 100% A
    p0[0] = 1

    T = int(tf/(tau/tau0))

    p  = msm.solve_FKE(p0, T)



#Plot time evolution of trimers and capsids

    plt.plot(np.array(range(T+1))*(tau/tau0),p[:,MM_dict['ndimer=90_nCD=30']],color=colors[cnt],label=r'MSM, $\tau/\tau_0=%.0f$' % (tau/tau0))
    #plt.plot(np.array(range(T+1)),p[:,MM_dict['ndimer=120_nCD=60']],color=colors[cnt],label=r'MSM, $\tau/\tau_0=%.0f$' % (tau/tau0))
    cnt += 1
plt.plot(range(tf),mc_yield,color='black',linestyle='--',label='MC')
plt.xlim([0,tf])
plt.ylim([0,0.1])
#plt.xscale('log')
#plt.xlim([100,3*10**4])
plt.xlabel(r'MC sweeps/$10^4$')
plt.ylabel(r'probability')
plt.legend(loc='lower right',fontsize=8,ncol=1)
plt.savefig('plots/T4_std_T3_yield_vs_time_mc_vs_msm_vary_tau.png')
