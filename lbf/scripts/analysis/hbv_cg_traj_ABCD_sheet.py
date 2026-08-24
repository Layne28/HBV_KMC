'''
Get a coarse-grained trajectory of HBV assembly.
States are defined via n_AB and n_CD_Hex
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
        myfile = myfolder + '/seed=%d/prod/energy_trunc.dat' % (i+1)
        if os.path.exists(myfile):
            #print(myfile)
            if os.path.getsize(myfile) > 0:
                print(myfile)
                process_trajectory(myfile)

def process_trajectory(trajectory_file):
    #Define intervals to use for defining states
    n_AB_interv = 1
    n_CD_Hex_interv = 1
    n_CD_T4_interv = 1

    #define max values of n_AB and n_CD (hex)
    n_AB_max=240
    n_CD_Hex_max=10#240
    n_CD_T4_max=240

    n_AB_states = math.ceil(n_AB_max/n_AB_interv)
    n_CD_Hex_states = math.ceil(n_CD_Hex_max/n_CD_Hex_interv)
    n_CD_T4_states = math.ceil(n_CD_T4_max/n_CD_T4_interv)

    bins_AB = []
    bins_CD_Hex = []
    bins_CD_T4 = []

    for i in range(n_AB_states):
        if n_AB_interv==1:
            n_AB_string = '%d' % i
        else:
            n_AB_string = '%d-%d' % (i*n_AB_interv, i*n_AB_interv+n_AB_interv-1)
        bins_AB.append(n_AB_string)

    for i in range(n_CD_Hex_states):
        if n_CD_Hex_interv==1:
            n_CD_Hex_string = '%d' % i
        else:
            n_CD_Hex_string = '%d-%d' % (i*n_CD_Hex_interv, i*n_CD_Hex_interv+n_CD_Hex_interv-1)
        bins_CD_Hex.append(n_CD_Hex_string)

    for i in range(n_CD_T4_states):
        if n_CD_T4_interv==1:
            n_CD_T4_string = '%d' % i
        else:
            n_CD_T4_string = '%d-%d' % (i*n_CD_T4_interv, i*n_CD_T4_interv+n_CD_T4_interv-1)
        bins_CD_T4.append(n_CD_T4_string)

    #Make a list of possible states
    state_list = []
    for i in range(n_AB_states):
        for j in range(n_CD_Hex_states):
            for k in range(n_CD_T4_states):
                if n_AB_interv==1:
                    n_AB_string = '%d' % i
                else:
                    n_AB_string = '%d-%d' % (i*n_AB_interv, i*n_AB_interv+n_AB_interv-1)
                if n_CD_Hex_interv==1:
                    n_CD_Hex_string = '%d' % j
                else:
                    n_CD_Hex_string = '%d-%d' % (j*n_CD_Hex_interv, j*n_CD_Hex_interv+n_CD_Hex_interv-1)
                if n_CD_T4_interv==1:
                    n_CD_T4_string = '%d' % k
                else:
                    n_CD_T4_string = '%d-%d' % (k*n_CD_T4_interv, k*n_CD_T4_interv+n_CD_T4_interv-1)
                state_list.append("nAB=%s_nCDHex=%s_nCDT4=%s" % (n_AB_string, n_CD_Hex_string, n_CD_T4_string))
    #state_list.append('sheet')

    #print(state_list)

    #Read in trajectory file to coarse-grain
    myfile = trajectory_file#sys.argv[1] #e.g. energy.dat
    myfolder = '/'.join(myfile.split('/')[:-1])
    print(myfolder)

    data_time = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols=('sweep'))
    data_nAB = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols=('NAB'))
    data_nCD =  np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols=('NCD_Hex','NCD_other','NCD_T4','NCD_T3'))

    #Convert to numpy arrays
    num_time = np.array(data_time['sweep'])
    num_nAB = np.array(data_nAB['NAB'])
    num_nCD_Hex = np.array(data_nCD['NCD_Hex'])
    num_nCD_T4 = np.array(data_nCD['NCD_T4'])
    #num_nCD = np.array([data_nCD['NCD_Hex'], data_nCD['NCD_other'], data_nCD['NCD_T4'], data_nCD['NCD_T3']]).T
    #num_nCD = np.sum(num_nCD, axis=1)

    #Remove rows with duplicate times
    data_tot = np.c_[num_time, num_nAB, num_nCD_Hex, num_nCD_T4]
    data_tot = np.unique(data_tot, axis=0)
    num_time = data_tot[:,0]
    num_nAB = data_tot[:,1]
    num_nCD_Hex = data_tot[:,2]
    num_nCD_T4 = data_tot[:,3]

    #Convert to states
    AB_states = []
    CD_Hex_states = []
    CD_T4_states = []
    is_malformed = []
    for i in range(num_nAB.shape[0]):
        malformed_flag = 0
        nAB = num_nAB[i]
        index = int(nAB//n_AB_interv)
        if index>=len(bins_AB):
            index = len(bins_AB)-1
        AB_states.append(bins_AB[index])

        nCD_Hex = num_nCD_Hex[i]
        index = int(nCD_Hex//n_CD_Hex_interv)
        if index>=len(bins_CD_Hex):
            index = len(bins_CD_Hex)-1
        if nCD_Hex>=n_CD_Hex_max:
            malformed_flag = 1
        CD_Hex_states.append(bins_CD_Hex[index])

        nCD_T4 = num_nCD_T4[i]
        index = int(nCD_T4//n_CD_T4_interv)
        if index>=len(bins_CD_T4):
            index = len(bins_CD_T4)-1
        CD_T4_states.append(bins_CD_T4[index])

        is_malformed.append(malformed_flag)
    
    state_traj = []
    for i in range(len(AB_states)):
        if is_malformed[i] == 1:
            state_traj.append("sheet")
        else:
            state_traj.append("nAB=%s_nCDHex=%s_nCDT4=%s" % (AB_states[i], CD_Hex_states[i], CD_T4_states[i]))
        #state_traj.append("nAB=%s_nCDHex=%s_nCDT4=%s" % (AB_states[i], CD_Hex_states[i], CD_T4_states[i]))

    #print(state_traj)

    with open(myfolder + '/cg_traj_ABCD_sheet.txt', 'w') as f:
        f.write('# time macrostate\n')
        for i in range(len(state_traj)):
            csv.writer(f, delimiter= " ", lineterminator="\n").writerow([num_time[i], state_traj[i]])

if __name__=="__main__":
    main()
