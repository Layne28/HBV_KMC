import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import matplotlib.image as image
import sys
import os

myfile = sys.argv[1] #energy.dat
column2 = 'Nd'

fig = plt.figure()

data = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True, usecols=('sweep','NE',column2)) 
data = np.array(data)
plt.plot(data['sweep'], data['NE'],color='blue',label='dimers')
plt.plot(data['sweep'], data[column2],color='red',label='CAMs',linewidth=0.8)
plt.xlabel('MC sweeps')
plt.ylabel('no. of dimers/CAMs')
#plt.ylim([0,150])
#plt.xlim([0,10*10**6])
plt.legend()

# #insert images
# im1 = image.imread('low_salt_trimer.tga')
# imagebox = OffsetImage(im1, zoom=0.04)
# ab = AnnotationBbox(imagebox, xy=(1*10**6,20),xybox=(0.7*10**6,25),frameon=False, bboxprops=dict(edgecolor='black'),arrowprops=dict(arrowstyle='-',color='black'))
# ax = plt.gca()
# ax.add_artist(ab)

# im2 = image.imread('low_salt_9mer.tga')
# imagebox = OffsetImage(im2, zoom=0.04)
# ab = AnnotationBbox(imagebox, xy=(2*10**6,30),xybox=(2.05*10**6,35),frameon=False, bboxprops=dict(edgecolor='black'),arrowprops=dict(arrowstyle='-',color='black'))
# ax = plt.gca()
# ax.add_artist(ab)

# im3 = image.imread('low_salt_15mer.tga')
# imagebox = OffsetImage(im3, zoom=0.06)
# ab = AnnotationBbox(imagebox, xy=(3.5*10**6,50),xybox=(3.5*10**6,55),frameon=False, bboxprops=dict(edgecolor='black'),arrowprops=dict(arrowstyle='-',color='black'))
# ax = plt.gca()
# ax.add_artist(ab)

# #Add markers for the 'CD-CD-CD' triangles
# plt.scatter(0.7*10**6,22.5,marker='o',color='purple',s=15,zorder=10)

# plt.scatter(1.88*10**6,29,marker='o',color='purple',s=15,zorder=10)
# plt.scatter(2.25*10**6,29,marker='o',color='purple',s=15,zorder=10)
# plt.scatter(2.08*10**6,40,marker='o',color='purple',s=15,zorder=10)

# plt.scatter(3.2*10**6,42,marker='o',color='purple',s=15,zorder=10)
# plt.scatter(3.6*10**6,42,marker='o',color='purple',s=15,zorder=10)
# plt.scatter(3.4*10**6,55,marker='o',color='purple',s=15,zorder=10)
# plt.scatter(3.8*10**6,55,marker='o',color='purple',s=15,zorder=10)
# plt.scatter(3.6*10**6,67,marker='o',color='purple',s=15,zorder=10)

plt.savefig('plots/ne_other_vs_time.png')
#plt.show()