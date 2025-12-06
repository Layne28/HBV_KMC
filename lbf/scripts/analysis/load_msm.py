import pickle
import sys

myfile = sys.argv[1] #e.g. msm_tau=1.000000.pkl

with open(myfile, 'rb') as f:
    msm = pickle.load(f)

print(msm.__dict__)
print(msm.compute_committor([0],[20]))
MM = msm.__dict__['_MSM__macrostate_map']
print(MM.__dict__)