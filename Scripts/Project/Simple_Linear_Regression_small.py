from sklearn.linear_model import LinearRegression
from sklearn.model_selection import RepeatedKFold
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error as mae

import random
import numpy as np

import matplotlib.pyplot as plt

random.seed(73)

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



def make_Linear_Regression(kf, normalized_matrix, Q3_percentage):
    """makes a 10 fold cv with a linear regression model"""
    
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
        #predictions2 = reg.predict(normalized_matrix[train_index]).ravel() ##

        #plt.scatter(Q3_percentage[test_index],list(predictions))
        #plt.scatter(Q3_percentage[train_index],list(predictions2))
        #plt.show()
        

        predictions = predictions.ravel()
        predictions_.append(predictions)
        
        measurements = Q3_percentage[test_index].ravel()
        measurements_.append(measurements)

        pearson_corr = pearsonr(measurements, predictions)
        rmse = mae(measurements, predictions, squared=False)

        
        print("R2: ", r2)
        print("Pearson Corr.: ", pearson_corr[0])
        print("RMSE: ", rmse)
        print("--------------------------------")
        print()

        r2_.append(r2)
        pearson_corr_.append(pearson_corr[0])
        rmse_.append(rmse)


    print(sum(pearson_corr_)/len(pearson_corr_))
    return np.concatenate(predictions_).ravel(), np.concatenate(measurements_).ravel()


def get_data_subset(counts_peaks_matrixes, Q3_percentage):
    """gets subset of data"""
    random_samples = random.sample(range(len(counts_peaks_matrixes)), 100)
    subset_counts_peaks_matrixes = [counts_peaks_matrixes[num] for num in random_samples]
    subset_Q3_percentage = [Q3_percentage[num] for num in random_samples]

    helix_percentage = [[sec[2]] for sec in subset_Q3_percentage]
    sheet_percentage = [[sec[1]] for sec in subset_Q3_percentage]
    coil_percentage = [[sec[0]] for sec in subset_Q3_percentage]

    return subset_counts_peaks_matrixes, helix_percentage, sheet_percentage, coil_percentage

def prepare_matrix(subset_counts_peaks_matrixes):
    """makes an array and normalizes it"""
    reshaped_matrix = reshape_to_1D(subset_counts_peaks_matrixes)
    normalized_matrix = normalization(reshaped_matrix)
    return normalized_matrix
    
if __name__=="__main__":
    with open("peak_matrixes_10x10_522.npy", "rb") as infile:
        counts_peaks_matrixes = np.load(infile)
        print(len(counts_peaks_matrixes))

    with open("Q3_percentage_522.npy", "rb") as infile:
        Q3_percentage = np.load(infile)

    print("SIMPLE LINEAR REGRESSION SMALL")
    for i in range(1):
        print(i, ". run ...")
        subset_counts_peaks_matrixes, helix_percentage, sheet_percentage, coil_percentage = get_data_subset(counts_peaks_matrixes, Q3_percentage)
        normalized_matrix = prepare_matrix(subset_counts_peaks_matrixes)

   
        rkf = RepeatedKFold(n_splits=4, n_repeats=4)
        #make_Linear_Regression(rkf, np.array(normalized_matrix), np.array(helix_percentage))


        #reg = LinearRegression().fit(normalized_matrix, np.array(helix_percentage))
        #predictions = reg.predict(normalized_matrix)
        #plt.scatter(helix_percentage, predictions)
        #plt.show()
        
