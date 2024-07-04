import sys

file = sys.argv[1] #expects 'trajlammps.dat'

out_folder = '/'.join((file.split('/'))[:-1])
print((file.split('/'))[:-1])
print(out_folder)

with open(file,'r') as f:
    a =f.readlines()
indices=[i for i, j in enumerate(a) if 'LAMMPS' in j]
indices.append(len(a))
for i in range(len(indices)-1):
    initial=indices[i]
    final=indices[i+1]
    with open(out_folder + '/frame%s.dat'%i,'w') as fi:
        for index in range(initial,final):
            fi.write(a[index])
print('%d files'%(len(indices)-1))
