'''
Truncate trajectories at some well-defined stopping points
-- NE>=250
-- Growth slows to some speed over some time interval (once assembly is big enough) 
'''

import numpy as np
import sys

myfile = sys.argv[1] #/path/to/energy.dat
base_folder = myfile.split('energy.dat')[0]
outfile = base_folder + 'energy_trunc.dat'

cut_index = -1
npace = 0
avgpace = 0
avgAddInterval=10000
lastNhe = 0
lastNheGrowth = 0
check_factor = 100 #1000 This sets the cutoff growth speed
dimer_add_rate = 10000

data = np.genfromtxt(myfile, dtype=float, delimiter=',', names=True)
#print(data.shape)
#print(data)
my_dtype = data.dtype
header = ','.join(data.dtype.names)

print('final size:', data['NE'][-1])
if np.any(data['NE'] >= 240):
    indices = np.argwhere(data['NE'] >= 240)
    indices = indices.flatten()
    print('final structure >=240')
    print(indices)
    cut_index = indices[0]
    print('final size:', data['NE'][cut_index])

for i in range(data['sweep'].shape[0]) :
    #measure speed during typical 'growth phase'
    if data['sweep'][i] % 10000 == 0 and data['NE'][i]>35 and data['NE'][i]<75:
        #print(data['sweep'][i])
        thispace = (data['NE'][i]*2-lastNhe)/10000
        #print('thispace:', thispace)
        avgpace = (npace*avgpace+thispace)/(npace+1.0)
        avgAddInterval = pow(10, (-1*np.floor(np.log10(avgpace))))
        npace += 1
        lastNhe = data['NE'][i]*2
        #print('avgpace:', avgpace)
        #print('avgAddInterval:', avgAddInterval)
    #Compute rate of dimer addition 
    #if data['sweep'][i]>(check_factor*avgAddInterval):
        #ndimer_added = data['NE'][i]-data['NE'][i-int(check_factor*avgAddInterval/200)]
        #print(data['sweep'][i], ndimer_added, int(check_factor*avgAddInterval))
    '''
    if(data['NE'][i]>75 and i%(int(check_factor*avgAddInterval/200))==0):#check_factor*avgAddInterval):
        print(avgAddInterval)
        if(np.abs(2*data['NE'][i]-lastNhe)<=4):
            print('Assembly is >75 edges and growth has slowed to <=2 edges per', avgAddInterval*check_factor, '. Stopping.')
            print('Compare current sweep to max sweep:', data['sweep'][i], data['sweep'][-1])
            cut_index = i
            print('real final size:', data['NE'][cut_index])
            #print('Compare current index to size of array:', i, data['sweep'].shape[0]-2)
            break
        lastNhe = data['NE'][i]*2
    '''
        
    if data['NE'][i]>75 and data['sweep'][i]%(check_factor*avgAddInterval)==0:
        #print('TEST:', data['NE'][i], lastNheGrowth)
        #print('sweep:', data['sweep'][i])
        #print('rate:', check_factor*avgAddInterval)
        #print(np.abs(data['NE'][i]))
        #print(lastNheGrowth/2)
        if (np.abs(data['NE'][i]*2-lastNheGrowth)<=4):
            print('Assembly is >75 edges and growth has slowed to <=2 edges per', check_factor*avgAddInterval, '. Stopping.')
            #print('Compare current sweep to max sweep:', data['sweep'][i], data['sweep'][-1])
            cut_index = i
            print('real final size:', data['NE'][cut_index])
            #print('Compare current index to size of array:', i, data['sweep'].shape[0]-2)
            break
        lastNheGrowth = data['NE'][i]*2
        lastNhe = data['NE'][i]*2

if cut_index != -1 and data['NE'][-1]==120 and data['NCD_T4'][-1]==60:
    print('WARNING: ENDING T4 TRAJECTORY PREMATURELY!')

if cut_index == -1 and not((data['NE'][-1]==120 and data['NCD_T4'][-1]==60) or (data['NE'][-1]==90 and data['NCD_T3'][-1]==30) or data['sweep'][-1]==200000000):
    print('ERROR: FAILED TO STOP STALLED TRAJECTORY!')
    exit()
print('cut_index:', cut_index)

if cut_index==-1:
    trunc_data = data
else:
    trunc_data = data[:(cut_index+1)]
#print(trunc_data['NE'])
#print(trunc_data)

my_fmt = ''
for i in range(len(data.dtype.names)):
    if i==2 or i==3:
        my_fmt += '%.05f,'
    else:
        my_fmt += '%d,'
my_fmt = my_fmt[:-1]
#print(my_fmt)

#Save truncated trajectory to new file

np.savetxt(outfile, trunc_data, header=header, delimiter=',', comments='', fmt=my_fmt)
