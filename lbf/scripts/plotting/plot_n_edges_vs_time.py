import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

myfile = sys.argv[1]

data = pd.read_csv(myfile)

print(data.keys())
print(data['sweep'])

fig = plt.figure()
plt.plot(data['sweep'], data['Nsurf'])
plt.show()