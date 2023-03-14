import numpy as np


def reshape_to_1D(matrix):
    """reshapes a matrix to a 100x1 array"""
    reshaped_matrix = np.reshape(matrix, (len(matrix),100))
    return reshaped_matrix


def normalization(reshaped_matrix):
    """normalizes the counts in the matrix to 1"""
    normalized_matrix = []
    for hsqc in reshaped_matrix:
        total_count = np.sum(hsqc)
        normalized_hsqc = []
        for grid in hsqc:
            normalized_grid = grid / total_count
            normalized_hsqc.append(normalized_grid)
        normalized_matrix.append(normalized_hsqc)

    return normalized_matrix


def seperate_secondary_structures(Q3_percentage):
    """[coil,sheet,helix],... -> [coil,...],[sheet,...],[helix,...]"""
        
    helix_percentage = [[sec[2]] for sec in Q3_percentage]
    sheet_percentage = [[sec[1]] for sec in Q3_percentage]
    coil_percentage = [[sec[0]] for sec in Q3_percentage]      

        
    return helix_percentage, sheet_percentage, coil_percentage
