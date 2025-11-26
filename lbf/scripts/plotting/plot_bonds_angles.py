'''Plot bond lengths and angles from the processed data.'''

import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_bond_lengths(data, header, output_dir):
    bond_indices = [i for i, h in enumerate(header) if h.startswith('l')]
    bond_lengths = data[:, bond_indices]
    
    #Plot the bond lengths versus sweeps
    plt.figure(figsize=(10, 6))
    for i, idx in enumerate(bond_indices):
        plt.plot(np.arange(bond_lengths.shape[0]), bond_lengths[:, i], label=header[idx])
    
    plt.xlabel('Sweep')
    plt.ylabel('Bond Length')
    plt.title('Bond Lengths Over Sweeps')
    plt.legend()
    plt.grid()
    plt.savefig(f'{output_dir}/bond_lengths.png')
    plt.show()
    
    #Plot a histogram of bond lengths
    plt.figure(figsize=(10, 6))
    colors = []
    for i, idx in enumerate(bond_indices):
        n, bins, patches = plt.hist(bond_lengths[:, i], bins=30, alpha=0.5, label=header[idx],density=True)
        colors.append(patches[0].get_facecolor())
    #Make gaussian distributions with mean and std dev
    #colored the same as the corresponding histogram
    means = np.mean(bond_lengths, axis=0)
    std_devs = np.std(bond_lengths, axis=0)
    for i, idx in enumerate(bond_indices):
        print(plt.gca().lines)
        x = np.linspace(means[i] - 3*std_devs[i], means[i] + 3*std_devs[i], 100)
        y = (1 / (std_devs[i] * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - means[i]) / std_devs[i])**2)
        plt.plot(x, y, label=f'{header[idx]} Gaussian', linestyle='--',color=colors[i])
    plt.xlabel('Bond Length')
    plt.ylabel('Frequency')
    plt.title('Histogram of Bond Lengths')
    plt.legend()
    plt.grid()
    plt.savefig(f'{output_dir}/bond_lengths_hist.png')
    plt.show()

def plot_angles(data, header, output_dir):
    angle_indices = [i for i, h in enumerate(header) if h.startswith('a')]
    angles = data[:, angle_indices]
    
    #Plot the binding angles vs sweeps
    plt.figure(figsize=(10, 6))
    for i, idx in enumerate(angle_indices):
        plt.plot(np.arange(angles.shape[0]), angles[:, i], label=header[idx])
    
    plt.xlabel('Sweep')
    plt.ylabel('Angle (radians)')
    plt.title('Angles Over Sweeps')
    plt.legend()
    plt.grid()
    plt.savefig(f'{output_dir}/angles.png')
    plt.show()
    
    #Plot a histogram of angles
    colors = []
    plt.figure(figsize=(10, 6))
    for i, idx in enumerate(angle_indices):
        n, bins, patches = plt.hist(angles[:, i], bins=30, alpha=0.5, label=header[idx],density=True)
        colors.append(patches[0].get_facecolor())
    #Make gaussian distributions with mean and std dev
    #colored the same as the corresponding histogram
    means = np.mean(angles, axis=0)
    std_devs = np.std(angles, axis=0)
    for i, idx in enumerate(angle_indices):
        x = np.linspace(means[i] - 3*std_devs[i], means[i] + 3*std_devs[i], 100)
        y = (1 / (std_devs[i] * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - means[i]) / std_devs[i])**2)
        plt.plot(x, y, label=f'{header[idx]} Gaussian', linestyle='--',color=colors[i])
    plt.xlabel('Angle (radians)')
    plt.ylabel('Frequency')
    plt.title('Histogram of Angles')
    plt.legend()
    plt.grid()
    plt.savefig(f'{output_dir}/angles_hist.png')
    plt.show()

def plot_dihedrals(data, header, output_dir):
    dihedral_indices = [i for i, h in enumerate(header) if h.startswith('d')]
    dihedrals = data[:, dihedral_indices]
    
    #Plot the dihedral angles vs sweeps
    plt.figure(figsize=(10, 6))
    for i, idx in enumerate(dihedral_indices):
        plt.plot(np.arange(dihedrals.shape[0]), dihedrals[:, i], label=header[idx])

    plt.xlabel('Sweep')
    plt.ylabel('Dihedral Angle (radians)')
    plt.title('Dihedral Angles Over Sweeps')
    plt.legend()
    plt.grid()
    plt.savefig(f'{output_dir}/dihedrals.png')
    plt.show()
    
    #Plot a histogram of dihedral angles
    colors = []
    plt.figure(figsize=(10, 6))
    for i, idx in enumerate(dihedral_indices):
        n, bins, patches = plt.hist(dihedrals[:, i], bins=30, alpha=0.5, label=header[idx],density=True)
        colors.append(patches[0].get_facecolor())
    #Make gaussian distributions with mean and std dev
    #colored the same as the corresponding histogram
    means = np.mean(dihedrals, axis=0)
    std_devs = np.std(dihedrals, axis=0)
    for i, idx in enumerate(dihedral_indices):
        x = np.linspace(means[i] - 3*std_devs[i], means[i] + 3*std_devs[i], 100)
        y = (1 / (std_devs[i] * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - means[i]) / std_devs[i])**2)
        plt.plot(x, y, label=f'{header[idx]} Gaussian', linestyle='--',color=colors[i])
    plt.xlabel('Dihedral Angle (radians)')
    plt.ylabel('Frequency')
    plt.title('Histogram of Dihedral Angles')
    plt.legend()
    plt.grid()
    plt.savefig(f'{output_dir}/dihedrals_hist.png')
    plt.show()
    
def main(input_file):
    data = np.loadtxt(input_file, delimiter=',', skiprows=1)
    header = np.genfromtxt(input_file, delimiter=',', max_rows=1, dtype=str)

    output_dir = input_file.rsplit('/', 1)[0]
    #Plot the data
    plot_bond_lengths(data, header,output_dir)
    plot_angles(data, header,output_dir)
    plot_dihedrals(data, header,output_dir)
    
    #Also write out two files containing the average and standard deviation
    #of each bond length, angle, and dihedral angle
    #Write this out to the same directory as the input file
    
    averages = np.mean(data, axis=0)
    std_devs = np.std(data, axis=0)
    with open(f'{output_dir}/averages.csv', 'w') as f:
        f.write(','.join(header) + '\n')
        f.write(','.join(map(str, averages)) + '\n')
    with open(f'{output_dir}/std_devs.csv', 'w') as f:
        f.write(','.join(header) + '\n')
        f.write(','.join(map(str, std_devs)) + '\n')
    
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_bonds_angles.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    main(input_file)
    print("Plots generated successfully.")