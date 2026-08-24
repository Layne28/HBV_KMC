from MSM import MSM
from state import MacrostateMap

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import csv
import os
#os.environ['OPENBLAS_NUM_THREADS'] = '1'
import pickle

num_pad = 100
target_state='nAB=60_nCD=0' #'sheet'

def get_counts(msm, MM, trajfile, maxsweep, do_pad, do_add_abs):

    #maxsweep is the max trajectory length
    #if trajectory file stops earlier than this
    #then pad the MSM count with transitions
    #going from last frame state to last frame state
        

    trajs = []
    traj0 = []
    times = []
    with open(trajfile, 'r') as f:
        my_reader = csv.reader(f, delimiter=' ')
        header = next(my_reader)
        #print(header)
        for row in my_reader:
            #print(row)
            times.append(float(row[0]))
            traj0.append(row[-1])
    trajs.append(traj0)

    #now that we have some data, let's add transition counts to the msm.
    #Need to "pad" trajectories shorter than the max trajectory length
    delta_t = times[1]-times[0] 
    maxframe = int(maxsweep/delta_t)
    print('maxframe:', maxframe)
    tau = msm.get_lag()
    lag = int(tau/delta_t)
    print('lag: ', lag)
    if tau<delta_t:
        print('WARNING: lag time is less than separation between frames. Setting lag time to separation between frames, %f.' % delta_t)
        lag = 1
    for traj in trajs:
        #Trajectories that end in T=4 or T=3 are stopped early.
        #To ensure that the MSM reflects the effectively irreversible
        #formation of these structures, "pad" these trajectories 
        if do_pad==1:
            if traj[-1] == 'nAB=60_nCD=0':
                for i in range(num_pad):
                    traj.append('nAB=60_nCD=0')
        #to prevent degenerate eigenvectors,
        #add initial state to end of trajectories that end
        #in absorbing state other than T4
        if do_add_abs==1 and traj[-1] != target_state:
            traj.append('nAB=2_nCDHex=0_nCDT4=0')
        for i in range(len(traj)-lag):
            a = MM.state_to_index(traj[i])
            b = MM.state_to_index(traj[i+lag])
            msm.add_count(a,b,1)
        sfinal = MM.state_to_index(traj[-1])
        #for i in range(len(traj)-lag,len(traj)):
        #    print(i)
        #    print(len(traj))
        #    print(traj[i])
        #    a = MM.state_to_index(traj[i])
        #    msm.add_count(a,sfinal)   
        #    print('added %s to %s count at step %d' % (traj[i], traj[-1], i))
             
    return



def main():
    #construct an MSM and perform calculations

    #cg_traj_file = sys.argv[1] #macrostate trajectory
    cg_traj_folder = sys.argv[1] #folder containing macrostate trajectories
    #ndimer_int = int(sys.argv[2])
    #nCD_int = int(sys.argv[3])
    #ndrug_int = int(sys.argv[4])
    #tau = float(sys.argv[5]) #lag time
    tau = float(sys.argv[2]) #lag time
    do_pad = int(sys.argv[3]) #pad T=3, T=4 trajectories
    do_add_abs = int(sys.argv[4]) #add absorbing state to end of non-T=4 trajectories
    target_state = sys.argv[5]

    #check if already exists
    out_file = cg_traj_folder + '/msm_ABCD_do_pad=%d_do_add_abs=%d_tau=%f.pkl' % (do_pad, do_add_abs, tau)
    if os.path.isfile(out_file):
        print('MSM file already exists. Exiting.')
        return

    
    #Create an empty macrostate map object
    MM = MacrostateMap()

    #Get states from trajectories
    state_list = []
    subfolders = [e for e in os.listdir(cg_traj_folder) if e.startswith('seed')]
    nseeds = len(subfolders)
    print('num seeds:', nseeds)
    #Need to "pad" trajectories shorter than the max trajectory length
    maxsweep = 0
    skip_seeds = []
    #if 'T4_std' in cg_traj_folder:
    #    skip_seeds = [41, 57]
    #if 'CAM_high_salt' in cg_traj_folder:
    #    skip_seeds = [75]
    #if 'CAM_mod_salt' in cg_traj_folder:
    #    skip_seeds = [87]
    #if 'CAM_low_salt' in cg_traj_folder:
    #    skip_seeds = [54, 8]
    for i in range(nseeds):
        cg_traj_file = cg_traj_folder + '/seed=%d/prod/cg_traj_ABCD.txt' % (i+1)
        print(cg_traj_file)
        if os.path.exists(cg_traj_file) and (i+1) not in skip_seeds:
            with open(cg_traj_file) as f:
                lines = f.readlines()
                lines = lines[1:]
            if float(lines[-1].split(' ')[0])>maxsweep:
                maxsweep = float(lines[-1].split(' ')[0])
            for line in lines:
                line = line.rstrip()
                state = line.split(' ')[-1]
                #print(state)
                if state not in state_list:
                    state_list.append(state)
    if 'CAM_low_salt' in cg_traj_folder:
        maxsweep = 2*10**6
    else:
        maxsweep = 2*10**8
    print(maxsweep)
    for state in state_list:
        MM.update_maps(state)

    print(state_list)

    #Next, you can construct an empty MSM with a given lag time over these states
    num_states = MM.get_num_states() #get the number of states in the state space
    msm = MSM(num_states, lag = tau)

    #Get counts from coarse-grained trajectories
    for i in range(nseeds):
        cg_traj_file = cg_traj_folder + '/seed=%d/prod/cg_traj_ABCD.txt' % (i+1)
        print(cg_traj_file)
        if os.path.exists(cg_traj_file) and (i+1) not in skip_seeds:
            get_counts(msm, MM, cg_traj_file, maxsweep, do_pad, do_add_abs)

    #the count matrix is built. Now we finalize it to construct the transition matrices
    msm.finalize_counts(MM)

    #lets print out the count and transition matrices. The transition matrix should be
    #quite close to the one defined above in sample_trajectory
    print('counts:')
    print(msm.get_count_matrix())
        
    #Save MSM to file
    pickle_file = cg_traj_folder + '/msm_ABCD_do_pad=%d_do_add_abs=%d_target_state=%s_tau=%f.pkl' % (do_pad, do_add_abs, target_state, msm.get_lag())
    with open(pickle_file, 'wb') as f:
        pickle.dump(msm, f)

    return



if __name__ == "__main__":
    main()
