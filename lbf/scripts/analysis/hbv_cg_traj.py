'''
Get a coarse-grained trajectory of HBV assembly.
States are defined via n_dimer and n_CD
'''

import math
import sys
import numpy as np
import csv
import os

def main():
    myfolder = sys.argv[1]
    subfolders = [e for e in os.listdir(myfolder) if e.startswith('seed')]
    nseeds = len(subfolders)
    print(myfolder)
    for i in range(nseeds):
        #myfile = myfolder + '/seed-%d/energy.dat' % (i+1)
        myfile = myfolder + '/seed=%d/prod/energy.dat' % (i+1)
        if os.path.exists(myfile):
            #print(myfile)
            if os.path.getsize(myfile) > 0:
                print(myfile)
                process_trajectory(myfile)

def process_trajectory(trajectory_file):
    #Define intervals to use for defining states
    n_dimer_interv = 1
    n_CD_interv = 1

    #define max values of n_dimer and n_CD
    n_dimer_max=130
    n_CD_max=130

    n_dimer_states = math.ceil(n_dimer_max/n_dimer_interv)
    n_CD_states = math.ceil(n_CD_max/n_CD_interv)

    bins_dimer = []
    bins_CD = []

    for i in range(n_dimer_states):
        if n_dimer_interv==1:
            n_dimer_string = '%d' % i
        else:
            n_dimer_string = '%d-%d' % (i*n_dimer_interv, i*n_dimer_interv+n_dimer_interv-1)
        bins_dimer.append(n_dimer_string)

    for i in range(n_CD_states):
        if n_CD_interv==1:
            n_CD_string = '%d' % i
        else:
            n_CD_string = '%d-%d' % (i*n_CD_interv, i*n_CD_interv+n_CD_interv-1)
        bins_CD.append(n_CD_string)

    #Make a list of possible states
    state_list = []
    for i in range(n_dimer_states):
        for j in range(n_CD_states):
            if n_dimer_interv==1:
                n_dimer_string = '%d' % i
            else:
                n_dimer_string = '%d-%d' % (i*n_dimer_interv, i*n_dimer_interv+n_dimer_interv-1)
            if n_CD_interv==1:
                n_CD_string = '%d' % j
            else:
                n_CD_string = '%d-%d' % (j*n_CD_interv, j*n_CD_interv+n_CD_interv-1)
            state_list.append("ndimer=%s_nCD=%s" % (n_dimer_string, n_CD_string))

    #print(state_list)

    #Read in trajectory file to coarse-grain
    myfile = trajectory_file#sys.argv[1] #e.g. energy.dat
    myfolder = '/'.join(myfile.split('/')[:-1])
    print(myfolder)

    data_time = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols=('sweep'))
    data_ndimer = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols=('NE'))
    data_nCD =  np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols=('NCD_Hex','NCD_other','NCD_T4','NCD_T3'))

    #Convert to numpy arrays
    num_time = np.array(data_time['sweep'])
    num_ndimer = np.array(data_ndimer['NE'])
    num_nCD = np.array(data_nCD['NCD_T4'])
    #num_nCD = np.array([data_nCD['NCD_Hex'], data_nCD['NCD_other'], data_nCD['NCD_T4'], data_nCD['NCD_T3']]).T
    #num_nCD = np.sum(num_nCD, axis=1)

    #Remove rows with duplicate times
    data_tot = np.c_[num_time, num_ndimer, num_nCD]
    data_tot = np.unique(data_tot, axis=0)
    num_time = data_tot[:,0]
    num_ndimer = data_tot[:,1]
    num_nCD = data_tot[:,2]

    #Convert to states
    dimer_states = []
    CD_states = []
    for i in range(num_ndimer.shape[0]):
        ndimer = num_ndimer[i]
        index = int(ndimer//n_dimer_interv)
        if index>=len(bins_dimer):
            index = len(bins_dimer)-1
        dimer_states.append(bins_dimer[index])

        nCD = num_nCD[i]
        index = int(nCD//n_CD_interv)
        if index>=len(bins_CD):
            index = len(bins_CD)-1
        CD_states.append(bins_CD[index])
    
    state_traj = []
    for i in range(len(dimer_states)):
        state_traj.append("ndimer=%s_nCD=%s" % (dimer_states[i], CD_states[i]))

    #print(state_traj)

    with open(myfolder + '/cg_traj.txt', 'w') as f:
        f.write('# time macrostate\n')
        for i in range(len(state_traj)):
            csv.writer(f, delimiter= " ", lineterminator="\n").writerow([num_time[i], state_traj[i]])

if __name__=="__main__":
    main()
