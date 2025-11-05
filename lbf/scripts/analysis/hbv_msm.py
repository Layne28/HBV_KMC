from MSM import MSM
from state import MacrostateMap

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import csv
import os
import pickle

def get_counts(msm, MM, trajfile):

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
    delta_t = times[1]-times[0] 
    tau = msm.get_lag()
    lag = int(tau/delta_t)
    if tau<delta_t:
        print('WARNING: lag time is less than separation between frames. Setting lag time to separation between frames, %f.' % delta_t)
        lag = 1
    for traj in trajs:
        for i in range(len(traj)-lag):
            a = MM.state_to_index(traj[i])
            b = MM.state_to_index(traj[i+lag])
            msm.add_count(a,b,1)
    
    return



def main():
    #construct an MSM and perform calculations

    #cg_traj_file = sys.argv[1] #macrostate trajectory
    cg_traj_folder = sys.argv[1] #folder containing macrostate trajectories
    tau = float(sys.argv[2]) #lag time

    #Create an empty macrostate map object
    MM = MacrostateMap()

    #Get states from trajectory
    state_list = []
    subfolders = [e for e in os.listdir(cg_traj_folder) if e.startswith('seed')]
    nseeds = len(subfolders)
    for i in range(nseeds):
        #cg_traj_file = cg_traj_folder + '/seed=%d/prod/cg_traj.txt' % (i+1)
        cg_traj_file = cg_traj_folder + '/seed-%d/cg_traj.txt' % (i+1)
        print(cg_traj_file)
        if os.path.exists(cg_traj_file):
            with open(cg_traj_file) as f:
                lines = f.readlines()
                lines = lines[1:]
            for line in lines:
                line = line.rstrip()
                state = line.split(' ')[-1]
                #print(state)
                if state not in state_list:
                    state_list.append(state)

    for state in state_list:
        MM.update_maps(state)

    print(state_list)

    #Next, you can construct an empty MSM with a given lag time over these states
    num_states = MM.get_num_states() #get the number of states in the state space
    msm = MSM(num_states, lag = tau)

    #Get counts from coarse-grained trajectories
    for i in range(nseeds):
        #cg_traj_file = cg_traj_folder + '/seed=%d/prod/cg_traj.txt' % (i+1)
        cg_traj_file = cg_traj_folder + '/seed-%d/cg_traj.txt' % (i+1)
        print(cg_traj_file)
        if os.path.exists(cg_traj_file):
            get_counts(msm, MM, cg_traj_file)

    #the count matrix is built. Now we finalize it to construct the transition matrices
    msm.finalize_counts(MM)

    #lets print out the count and transition matrices. The transition matrix should be
    #quite close to the one defined above in sample_trajectory
    print('counts:')
    print(msm.get_count_matrix())
    #print(msm.get_transition_matrix())

    print('committor:')
    print(MM.state_to_index('ndimer=120_nCD=60'))

    print(msm.compute_committor([MM.state_to_index('ndimer=3_nCD=1'),MM.state_to_index('ndimer=3_nCD=3')],[MM.state_to_index('ndimer=120_nCD=60')]))

    #Now we can use the transition matrix to compute things. Lets solve the forward 
    #kolmogorov equation to predict the yield of monomer and dimers if we start with all
    #monomers
    cmap = mpl.colormaps['plasma']
    colors = cmap(np.linspace(0,1,msm._MSM__num_states))
    p0 = [0]*msm._MSM__num_states#[1,0] #initial distribution. 100% A
    p0[0] = 1
    T  = 20000  #final time (in lags) 
    p  = msm.solve_FKE(p0, T)
    print(p.shape)
    #plt.plot(range(T+1), p, linewidth=2, color=colors)
    for i in range(p.shape[1]):
        plt.plot(range(T+1), p[:,i], linewidth=2, color=colors[i])
    plt.xlabel("Lag Times")
    plt.ylabel("Yield")
    plt.legend(np.arange(msm._MSM__num_states).tolist())
    #plt.legend(["A", "B"])
    plt.show()

    # #if we start with all dimers instead, we just change p0
    # p0 = [0, 1] #initial distribution, 100% dimer
    # p  = msm.solve_FKE(p0, T)
    # plt.plot(range(T+1), p, linewidth=2)
    # plt.xlabel("Lag Times")
    # plt.ylabel("Yield")
    # plt.legend(["A", "B"])
    # plt.show()

    # #if we start with a 50/50 split of monomer and dimer
    # p0 = [0.5, 0.5] #initial distribution, 100% dimer
    # p  = msm.solve_FKE(p0, T)
    # plt.plot(range(T+1), p, linewidth=2)
    # plt.xlabel("Lag Times")
    # plt.ylabel("Yield")
    # plt.legend(["A", "B"])
    # plt.show()
    
    #Save MSM to file
    pickle_file = cg_traj_folder + '/msm_tau=%f.pkl' % msm.get_lag()
    with open(pickle_file, 'wb') as f:
        pickle.dump(msm, f)

    return



if __name__ == "__main__":
    main()