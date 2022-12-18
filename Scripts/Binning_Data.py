import pandas as pd
import numpy as np
import pickle


def binning(shifts, binsize, shift_min, num_1D_grid):
    """gives the indexes for certain bins by float dividing the shift value by the binsize"""
    bin_indexes = []
    for shifts_one_protein in shifts:
        bin_indexes_one_protein = []
        
        for ppm_value in shifts_one_protein:
            bin_index = int((ppm_value - shift_min) // binsize)
            if bin_index > (num_1D_grid - 1):
                #print(ppm_value, shift_min, binsize, "help me")
                bin_index = num_1D_grid - 1
            elif bin_index < 0:
                #print(ppm_value, shift_min, binsize, "what the fuck")
                bin_index = 0
                
            bin_indexes_one_protein.append(bin_index)

        bin_indexes.append(bin_indexes_one_protein)

    return bin_indexes


def get_shifts(chemical_shifts, atomtype):
    """returns a list of with all shifts of a certain atomtype. X stands for H and Y for N"""
    shifts = []
    
    for NMR_data in chemical_shifts.values():
        shifts_one_protein = NMR_data[f"{atomtype}_shift"].to_numpy()
        shifts.append(shifts_one_protein)

    return shifts
        

def generate_count_peaks_matrix(binned_H_shifts, binned_N_shifts, num_1D_grid):
    """counts peaks in a certain bin in the defined 2D grid"""
    grid_2D = (num_1D_grid, num_1D_grid)
    count_peaks_matrixes = []
    
    for index in range(len(binned_H_shifts)):
        
        count_peaks_matrix = np.zeros(grid_2D)
        #print(count_peaks_matrix)
    
        for H_shift_bin, N_shift_bin in zip(binned_H_shifts[index], binned_N_shifts[index]):
            count_peaks_matrix[N_shift_bin, H_shift_bin] += 1

        count_peaks_matrixes.append(count_peaks_matrix.copy())

    return count_peaks_matrixes


if __name__=='__main__':
    # you need also to read a file of the bmrb ids that were filtered out
    # and then take only those for the binning
    H_shift_max = 15
    H_shift_min = 5

    N_shift_max = 150
    N_shift_min = 80

    #num_1D_grid = int(input("Number of 1D grids: "))
    num_1D_grid = 5

    H_binsize = (H_shift_max - H_shift_min) / num_1D_grid
    N_binsize = (N_shift_max - N_shift_min) / num_1D_grid
    print(H_binsize)
    print(N_binsize)
    print("--------------")
    
    with open("chemical_shifts.pkl", "rb") as infile:
        chemical_shifts = pickle.load(infile)

    N_shifts = get_shifts(chemical_shifts, atomtype="Y")
    H_shifts = get_shifts(chemical_shifts, atomtype="X")

    #print(H_shifts[0])
    #print(N_shifts[0])

    binned_H_shifts = binning(H_shifts, H_binsize, H_shift_min, num_1D_grid)
    binned_N_shifts = binning(N_shifts, N_binsize, N_shift_min, num_1D_grid)

    print(binned_H_shifts[0])
    print(binned_N_shifts[0])

    count_peaks_matrixes = generate_count_peaks_matrix(binned_H_shifts, binned_N_shifts, num_1D_grid)
    print(count_peaks_matrixes[0])
