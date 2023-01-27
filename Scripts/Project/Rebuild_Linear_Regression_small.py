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

random.seed(73)

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
        

        reg = LinearRegression().fit(matrix, structure_percentage)
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

        #print(i)
        #print(len(train_index))
        #print(len(test_index))
        reg = LinearRegression().fit(normalized_matrix[train_index], Q3_percentage[train_index])
        r2 = reg.score(normalized_matrix[train_index], Q3_percentage[train_index])
        predictions = reg.predict(normalized_matrix[test_index])

        predictions = predictions.ravel()
        predictions_.append(predictions)
        
        measurements = Q3_percentage[test_index].ravel()
        measurements_.append(measurements)

        pearson_corr = pearsonr(measurements, predictions)
        rmse = mae(measurements, predictions, squared=False)

        
        #print("R2: ", r2)
        #print("Pearson Corr.: ", pearson_corr[0])
        #print("RMSE: ", rmse)
        #print("--------------------------------")
        #print()
 
    
    return np.concatenate(predictions_).ravel(), np.concatenate(measurements_).ravel()# currently only last predictions, if you make it as a list you might need to change the np.zeros number


                                        
def external_correction(preds, pearson_correlation):
    """correct"""
    
    external_preds = [list(), list(), list()]

    for index, (preds_X, pearson_X) in enumerate(zip(preds, pearson_correlation)):
        for P in preds_X:
            P_ext = P + ( P * (sum(pearson_X) / (sum(pearson_correlation[0]) + sum(pearson_correlation[1]) + sum(pearson_correlation[2]))  ) )
            #P_ext = P + ( P * ( (sum(pearson_correlation[0]) + sum(pearson_correlation[1]) + sum(pearson_correlation[2])) / sum(pearson_X)  ) )
            external_preds[index].append(P_ext)

    return external_preds


def internal_correction(preds, slope, intercept):
    """internal correction"""
    internal_preds = []
    for P in preds:
        P_int = (P - intercept) / slope
        internal_preds.append(P_int)

    return internal_preds


if __name__=="__main__":
    with open("peak_matrixes_10x10_522.npy", "rb") as infile:
        counts_peaks_matrixes = np.load(infile)
        print(len(counts_peaks_matrixes))

    with open("Q3_percentage_522.npy", "rb") as infile:
        Q3_percentage = np.load(infile)

    print("SIMPLE LINEAR REGRESSION SMALL")
    rkf = RepeatedKFold(n_splits=4, n_repeats=4)
    
    for i in range(1):
        print(i, ". run ...")
        subset_counts_peaks_matrixes, helix_percentage, sheet_percentage, coil_percentage = get_data_subset(counts_peaks_matrixes, Q3_percentage)
        structure_percentage = [helix_percentage, sheet_percentage, coil_percentage]
        normalized_matrix = prepare_matrix(subset_counts_peaks_matrixes)


        #
        grid_correlation_03_ = set()
        correlation_value_03_ = []
        
        for percentage in structure_percentage:
            bin_by_matrix_position = separate_bins(normalized_matrix)
            grid_correlation_03, correlation_value_03 = calculate_pearson_coefficient(bin_by_matrix_position, percentage)
            grid_correlation_03_ |= set(grid_correlation_03)
            correlation_value_03_.append(grid_correlation_03)

            
            
        preds_ = []
        measurements_ = []
        for percentage in structure_percentage:
            high_correlation_matrix_ = build_high_correlation_matrix(list(grid_correlation_03_), normalized_matrix)
            preds, measurements = linear_model(rkf, np.array(high_correlation_matrix_), np.array(percentage))

            preds_.append(preds)
            measurements_.append(measurements)

        
          
        external_preds = external_correction(preds_, correlation_value_03_)
        ##


        # pearson tests
        sum_pred = np.zeros(len(preds_[0]))
        for preds, measurements, ext_preds in zip(preds_, measurements_, external_preds):

            print(pearsonr(measurements, preds))
            print(mae(measurements, preds, squared=False))
            print()
            print(pearsonr(measurements, ext_preds))
            print(mae(measurements, ext_preds, squared=False))
            res = linregress(measurements, ext_preds)
            print("--------------------------------")

            #print(preds)
            #print(ext_preds)
            #print(measurements)
            #exit(0)


            internal_preds = internal_correction(ext_preds, res[0], res[1])

            print(pearsonr(measurements, internal_preds))
            print(mae(measurements, internal_preds, squared=False))
            print("###############################")
            print()
            print()
            sum_pred += np.array(preds)

            print()
            print()

        print(sum_pred)
            # things to do: now only last predictions, just writ them in a list
            # is it okay to draw a linear regression from all data --> corresponding to the paper yes
            # make the range stuff again
            # make plots

            
            
        
