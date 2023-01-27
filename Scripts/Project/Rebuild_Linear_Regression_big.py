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
from Rebuild_Linear_Regression_small import external_correction
from Rebuild_Linear_Regression_small import internal_correction

import scipy.stats
import warnings

warnings.filterwarnings("ignore", category=scipy.stats.ConstantInputWarning)

random.seed(73) #why doesn't it work!!!


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
        helix_percentage, sheet_percentage, coil_percentage = get_data_subset(Q3_percentage)
        structure_percentage = [helix_percentage, sheet_percentage, coil_percentage]
        normalized_matrix = prepare_matrix(counts_peaks_matrixes)


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

            
            
        
