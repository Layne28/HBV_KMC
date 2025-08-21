import matplotlib.pyplot as plt
import numpy as np

kappa = 40.0
phi0 = 0.75

fig = plt.figure()
xvals = np.linspace(0,1.1,100)
Evals = kappa*(1-np.cos(xvals-phi0))
yvals = np.exp(-Evals)

plt.plot(xvals, yvals)
plt.ylim([0,1])
plt.savefig('plots/dihedral_example.png')
plt.show()