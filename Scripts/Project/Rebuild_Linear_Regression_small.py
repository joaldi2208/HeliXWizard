from sklearn.linear_model import LinearRegression
from sklearn.model_selection import RepeatedKFold
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error as mae

from scipy.stats import linregress

import random
import numpy as np

from Simple_Linear_Regression_small import reshape_to_1D
from Simple_Linear_Regression_small import normalization
from Simple_Linear_Regression_small import get_data_subset
from Simple_Linear_Regression_small import prepare_matrix

import scipy.stats
import warnings
warnings.filterwarnings("ignore", category=scipy.stats.ConstantInputWarning)

import matplotlib.pyplot as plt

def separate_bins(normalized_matrix):
    """separate bins for the same molecule and connect them to bins for the same matrix position"""
    bin_by_matrix_position = [[[value] for value in matrix] for matrix in zip(*normalized_matrix)]
    return bin_by_matrix_position


def calculate_pearson_coefficient(bin_by_matrix_position, structure_percentage):
    """calculates the pearson coefficient for a linear regression"""

    grid_correlation_03 = []
    correlation_value_03 = []
    
    pearson_correlation_ = []
    
    for matrix_grid, matrix in enumerate(bin_by_matrix_position):
        

        reg = LinearRegression(fit_intercept=False).fit(matrix, structure_percentage)
        preds = reg.predict(matrix).ravel()
        pearson_correlation = pearsonr(preds, structure_percentage)
        
        
        if pearson_correlation[0] >= 0.3 or pearson_correlation[0] <= -0.3:
            grid_correlation_03.append(matrix_grid)
            correlation_value_03.append([pearson_correlation[0]])

                                        
    return grid_correlation_03, correlation_value_03



def build_high_correlation_matrix(grid_correlation_03, normalized_matrixes):
    """reduces matrix to only the values with high correlation"""
    high_correlation_matrix_ = []
    
    for matrix in normalized_matrixes:

        high_correlation_matrix = []
        for grid in grid_correlation_03:
            high_correlation_matrix.append(matrix[grid])

        high_correlation_matrix_.append(high_correlation_matrix)

    return high_correlation_matrix_

                                        
def linear_model(kf, normalized_matrix, Q3_percentage):
    """makes a linear regression and calculates the pearson correlation"""

    r2_ = []
    pearson_corr_ = []
    rmse_ = []
    predictions_ = []
    measurements_ = []

    
    for i, (train_index, test_index) in enumerate(kf.split(Q3_percentage)):

        #print(len(train_index))
        #print(len(test_index))
        reg = LinearRegression(fit_intercept=False).fit(normalized_matrix[train_index], Q3_percentage[train_index])
        r2 = reg.score(normalized_matrix[train_index], Q3_percentage[train_index])
        r2_.append(r2)
        predictions = reg.predict(normalized_matrix[test_index])

        predictions = predictions.ravel()
        predictions_.append(predictions)
        
        measurements = Q3_percentage[test_index].ravel()
        measurements_.append(measurements)

        pearson_corr = pearsonr(measurements, predictions)[0]
        pearson_corr_.append(pearson_corr)
        rmse = mae(measurements, predictions, squared=False)
        rmse_.append(rmse)
        
        #print("R2: ", r2)
        #print("Pearson Corr.: ", pearson_corr[0])
        #print("RMSE: ", rmse)
        #print("--------------------------------")
        #print()
 
    
    return np.concatenate(predictions_).ravel(), np.concatenate(measurements_).ravel(), pearson_corr_, rmse_, r2_



def rebuild_dataset_small():
    with open("peak_matrixes_10x10_522.npy", "rb") as infile:
        counts_peaks_matrixes = np.load(infile)

    with open("Q3_percentage_522.npy", "rb") as infile:
        Q3_percentage = np.load(infile)

    subset_counts_peaks_matrixes, helix_percentage, sheet_percentage, coil_percentage = get_data_subset(counts_peaks_matrixes, Q3_percentage)
    structure_percentage = [helix_percentage, sheet_percentage, coil_percentage]
    normalized_matrix = prepare_matrix(subset_counts_peaks_matrixes)

    grid_correlation_03_ = set()
    correlation_value_03_ = []
        
    for percentage in structure_percentage:
        bin_by_matrix_position = separate_bins(normalized_matrix)
        grid_correlation_03, correlation_value_03 = calculate_pearson_coefficient(bin_by_matrix_position, percentage)
        grid_correlation_03_ |= set(grid_correlation_03)
        correlation_value_03_.append(grid_correlation_03)

            
    high_correlation_matrix_ = build_high_correlation_matrix(list(grid_correlation_03_), normalized_matrix)

    return helix_percentage, sheet_percentage, coil_percentage, high_correlation_matrix_

def rebuild_linear_regression_small(sec_structure, high_correlation_matrix_):
    rkf = RepeatedKFold(n_splits=4, n_repeats=4)
    predictions, measurements, pearson, rmse, r2 = linear_model(rkf, np.array(high_correlation_matrix_), np.array(sec_structure))
        
    return measurements, predictions, pearson, rmse, r2
            
if __name__=="__main__":
    helix_percentage, sheet_percentage, coil_percentage, high_correlation_matrix = rebuild_dataset_small()
    h_pearson_ = []
    s_pearson_ = []
    c_pearson_ = []

    h_rmse_ = []
    s_rmse_ = []
    c_rmse_ = []

    h_predictions_ = []
    s_predictions_ = []
    c_predictions_ = []
    
    sum_predictions_ = []
    for i in range(1):
        h_measurements, h_predictions, h_pearson, h_rmse, h_r2 = rebuild_linear_regression_small(helix_percentage, high_correlation_matrix)
        s_measurements, s_predictions, s_pearson, s_rmse, s_r2 = rebuild_linear_regression_small(sheet_percentage, high_correlation_matrix)
        c_measurements, c_predictions, c_pearson, c_rmse, c_r2 = rebuild_linear_regression_small(coil_percentage, high_correlation_matrix)

        sum_predictions = h_predictions + s_predictions + c_predictions
        
        sum_predictions_ += list(sum_predictions)
        h_pearson_ += list(h_pearson)
        s_pearson_ += list(s_pearson)
        c_pearson_ += list(c_pearson)

        h_rmse_ += list(h_rmse)
        s_rmse_ += list(s_rmse)
        c_rmse_ += list(c_rmse)

        h_predictions_ += list(h_predictions)
        s_predictions_ += list(s_predictions)
        c_predictions_ += list(c_predictions)

        
        
    print(max(sum_predictions_), sum(sum_predictions_)/len(sum_predictions_), np.median(sum_predictions_), min(sum_predictions_))
    all_predictions = np.array(h_predictions_) + np.array(s_predictions_) + np.array(c_predictions_)
    print()
    print(np.average(all_predictions), np.median(all_predictions))
    print("---")

    print(np.average(h_pearson), np.median(h_pearson))
    fig, ax = plt.subplots(figsize=(10,6))
    ax.boxplot([h_pearson_, s_pearson, c_pearson])
    plt.xticks([1, 2, 3], ["helix", "sheet", "coil"])
    plt.show()

    print(np.average(h_rmse), np.median(h_rmse))
    fig, ax = plt.subplots(figsize=(10,6))
    ax.boxplot([h_rmse_, s_rmse, c_rmse])
    plt.xticks([1, 2, 3], ["helix", "sheet", "coil"])
    plt.show()
