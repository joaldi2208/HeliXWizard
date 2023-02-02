from sklearn.linear_model import LinearRegression
from sklearn.model_selection import RepeatedKFold
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error as mae

from scipy.stats import linregress

import random
import numpy as np

from Simple_Linear_Regression_big import reshape_to_1D
from Simple_Linear_Regression_big import normalization
from Simple_Linear_Regression_big import get_data_subset
from Simple_Linear_Regression_big import prepare_matrix

from Rebuild_Linear_Regression_small import separate_bins
from Rebuild_Linear_Regression_small import calculate_pearson_coefficient
from Rebuild_Linear_Regression_small import build_high_correlation_matrix
from Rebuild_Linear_Regression_small import linear_model

import matplotlib.pyplot as plt

import scipy.stats
import warnings

warnings.filterwarnings("ignore", category=scipy.stats.ConstantInputWarning)


def rebuild_dataset_big():
    with open("peak_matrixes_10x10_4087.npy", "rb") as infile:
        counts_peaks_matrixes = np.load(infile)

    with open("Q3_percentage_4087.npy", "rb") as infile:
        Q3_percentage = np.load(infile)

    helix_percentage, sheet_percentage, coil_percentage = get_data_subset(Q3_percentage)
    structure_percentage = [helix_percentage, sheet_percentage, coil_percentage]
    normalized_matrix = prepare_matrix(counts_peaks_matrixes)

    grid_correlation_03_ = set()
    correlation_value_03_ = []
    
    for percentage in structure_percentage:
        bin_by_matrix_position = separate_bins(normalized_matrix)
        grid_correlation_03, correlation_value_03 = calculate_pearson_coefficient(bin_by_matrix_position, percentage)
        grid_correlation_03_ |= set(grid_correlation_03)
        correlation_value_03_.append(grid_correlation_03)

    high_correlation_matrix_ = build_high_correlation_matrix(list(grid_correlation_03_), normalized_matrix)

    return helix_percentage, sheet_percentage, coil_percentage, high_correlation_matrix_


def rebuild_linear_regression_big(sec_structure, high_correlation_matrix_):
    rkf = RepeatedKFold(n_splits=4, n_repeats=4)
    predictions, measurements, pearson, rmse, r2 = linear_model(rkf, np.array(high_correlation_matrix_), np.array(sec_structure))

    return measurements, predictions, pearson, rmse, r2
            
            
if __name__=="__main__":
    helix_percentage, sheet_percentage, coil_percentage, high_correlation_matrix = rebuild_dataset_big()
    
    h_pearson_ = []
    s_pearson_ = []
    c_pearson_ = []

    h_rmse_ = []
    s_rmse_ = []
    c_rmse_ = []
    
    sum_predictions_ = []
    for i in range(1):
        h_measurements, h_predictions, h_pearson, h_rmse, h_r2 = rebuild_linear_regression_big(helix_percentage, high_correlation_matrix)
        s_measurements, s_predictions, s_pearson, s_rmse, s_r2 = rebuild_linear_regression_big(sheet_percentage, high_correlation_matrix)
        c_measurements, c_predictions, c_pearson, c_rmse, c_r2 = rebuild_linear_regression_big(coil_percentage, high_correlation_matrix)

        sum_predictions = h_predictions + s_predictions + c_predictions
        sum_predictions_ += list(sum_predictions)
        h_pearson_ += list(h_pearson)
        s_pearson_ += list(s_pearson)
        c_pearson_ += list(c_pearson)

        h_rmse_ += list(h_rmse)
        s_rmse_ += list(s_rmse)
        c_rmse_ += list(c_rmse)
        
    print(max(sum_predictions_), sum(sum_predictions_)/len(sum_predictions_), np.median(sum_predictions_), min(sum_predictions_))


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
