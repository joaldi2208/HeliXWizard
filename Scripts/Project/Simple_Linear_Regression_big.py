from sklearn.model_selection import RepeatedKFold
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error as mae

import random
import numpy as np

from Simple_Linear_Regression_small import make_Linear_Regression
#random.seed(73)

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



def get_data_subset(Q3_percentage):
    """gets subset of data"""

    helix_percentage = [[sec[2]] for sec in Q3_percentage]
    sheet_percentage = [[sec[1]] for sec in Q3_percentage]
    coil_percentage = [[sec[0]] for sec in Q3_percentage]

    return helix_percentage, sheet_percentage, coil_percentage

def prepare_matrix(subset_counts_peaks_matrixes):
    """makes an array and normalizes it"""
    reshaped_matrix = reshape_to_1D(subset_counts_peaks_matrixes)
    normalized_matrix = normalization(reshaped_matrix)
    return normalized_matrix
    

def simple_dataset_big():
    with open("peak_matrixes_10x10_4087.npy", "rb") as infile:
        counts_peaks_matrixes = np.load(infile)

    with open("Q3_percentage_4087.npy", "rb") as infile:
        Q3_percentage = np.load(infile)
        
    helix_percentage, sheet_percentage, coil_percentage = get_data_subset(Q3_percentage)
    normalized_matrix = prepare_matrix(counts_peaks_matrixes)

    return helix_percentage, sheet_percentage, coil_percentage, normalized_matrix

def simple_linear_regression_big(sec_structure, normalized_matrix):
    rkf = RepeatedKFold(n_splits=4, n_repeats=4)
    predictions, measurements, pearson, rmse, r2 = make_Linear_Regression(rkf, np.array(normalized_matrix), np.array(sec_structure))

    return measurements, predictions, pearson, rmse, r2

if __name__=="__main__":
    helix_percentage, sheet_percentage, coil_percentage, normalized_matrix = simple_dataset_big()
    for i in range(10):
        measurements, predictions, pearson, rmse, r2 = simple_linear_regression_big(helix_percentage, normalized_matrix)
        measurements, predictions, pearson, rmse, r2 = simple_linear_regression_big(sheet_percentage, normalized_matrix)
        measurements, predictions, pearson, rmse, r2 = simple_linear_regression_big(coil_percentage, normalized_matrix)
  
