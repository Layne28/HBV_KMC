from MSM import MSM
from state import MacrostateMap

import numpy as np
import matplotlib.pyplot as plt
import sys
import csv

def get_counts(msm, MM, trajfile):

    trajs = []
    traj0 = []
    times = []
    with open(trajfile, 'r') as f:
        my_reader = csv.reader(f, delimiter=' ')
        header = next(my_reader)
        print(header)
        for row in my_reader:
            print(row)
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

    cg_traj_file = sys.argv[1] #macrostate trajectory
    tau = float(sys.argv[2]) #lag time

    #Create an empty macrostate map object
    MM = MacrostateMap()

    #List states for a double-well barrier-crossing transition. 
    MM.update_maps("A") #add state A to maps
    MM.update_maps("B") #add state B to maps

    #Next, you can construct an empty MSM with a given lag time over these states
    num_states = MM.get_num_states() #get the number of states in the state space
    msm = MSM(num_states, lag = tau)

    #Get counts from coarse-grained trajectory
    get_counts(msm, MM, cg_traj_file)

    #the count matrix is built. Now we finalize it to construct the transition matrices
    msm.finalize_counts(MM)

    #lets print out the count and transition matrices. The transition matrix should be
    #quite close to the one defined above in sample_trajectory
    print(msm.get_count_matrix())
    print(msm.get_transition_matrix())

    #Now we can use the transition matrix to compute things. Lets solve the forward 
    #kolmogorov equation to predict the yield of monomer and dimers if we start with all
    #monomers
    p0 = [1,0] #initial distribution. 100% A
    T  = 200  #final time (in lags) 
    p  = msm.solve_FKE(p0, T)
    plt.plot(range(T+1), p, linewidth=2)
    plt.xlabel("Lag Times")
    plt.ylabel("Yield")
    plt.legend(["A", "B"])
    plt.show()

    #if we start with all dimers instead, we just change p0
    p0 = [0, 1] #initial distribution, 100% dimer
    p  = msm.solve_FKE(p0, T)
    plt.plot(range(T+1), p, linewidth=2)
    plt.xlabel("Lag Times")
    plt.ylabel("Yield")
    plt.legend(["A", "B"])
    plt.show()

    #if we start with a 50/50 split of monomer and dimer
    p0 = [0.5, 0.5] #initial distribution, 100% dimer
    p  = msm.solve_FKE(p0, T)
    plt.plot(range(T+1), p, linewidth=2)
    plt.xlabel("Lag Times")
    plt.ylabel("Yield")
    plt.legend(["A", "B"])
    plt.show()

    return



if __name__ == "__main__":
    main()