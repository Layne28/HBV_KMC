import numpy as np
import sys


def process_data(input_file, output_file):
    '''Read in length, angle, and dihedral data from the input file, remove duplicate and non-changing columns, and write the processed data to the output file.'''
    
    #Load header labels into numpy array
    with open(input_file, 'r') as f:
        header = f.readline().strip().split(',')
        header = np.array(header)
    
    data = np.loadtxt(input_file, delimiter=',', skiprows=1)
    data = data[:,3:] # Remove the first three columns (sweep, seed, seconds)
    header = header[3:] # Adjust header to match data
    
    # Remove unchanging columns
    header = header[~np.all(data == data[0, :], axis=0)]
    data = data[:, ~np.all(data == data[0, :], axis=0)]
    
    # Remove duplicate columns
    unique_indices = np.unique(data, axis=1, return_index=True)[1]
    header = header[np.sort(unique_indices)]
    data = data[:, np.sort(unique_indices)]
    
    # Remove columns with nan
    valid_columns = ~np.isnan(data).any(axis=0)
    header = header[valid_columns]
    data = data[:, valid_columns]
    
    #Now remove the columns with header a8 or a10
    unwanted_headers = ['a8', 'a10']
    unwanted_indices = np.isin(header, unwanted_headers)
    header = header[~unwanted_indices]
    data = data[:, ~unwanted_indices]
    
    #Finally, replace the half-edge headings with edge headings
    #l3 becomes l2
    #l5 becomes l3
    #l7 becomes l4
    #l9 becomes l5
    #a3 becomes a2
    #a5 becomes a3
    #a7 becomes a5
    #a9 becomes a6
    #d3 becomes d1
    header = np.array([h.replace('l3', 'l2').replace('l5', 'l3').replace('l7', 'l4').replace('l9', 'l5')
                       .replace('a3', 'a2').replace('a5', 'a3').replace('a7', 'a5').replace('a9', 'a6')
                       .replace('d3', 'd1') for h in header])

    #Swap columns a3 and a4
    header[[7,8]] = header[[8,7]]
    data[:,[7,8]] = data[:,[8,7]]  

    print(header)
    print(data.shape)
    # Write the processed data to the output file
    header = ','.join(header)
    print(header)
    
    # Save the processed data to the output file
    np.savetxt(output_file, data, header=header, delimiter=',',comments='')
    
def main():
    if len(sys.argv) != 3:
        print("Usage: python postprocess_angles.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    process_data(input_file, output_file)
    print(f"Processed data saved to {output_file}")

if __name__ == "__main__":
    main()