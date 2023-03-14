from sklearn.linear_model import LinearRegression
import scipy.stats
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=scipy.stats.ConstantInputWarning)

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
        pearson_correlation = scipy.stats.pearsonr(preds, structure_percentage)

        
        if pearson_correlation[0] >= 0.3 or pearson_correlation[0] <= -0.3:
            grid_correlation_03.append(matrix_grid)
            correlation_value_03.append([pearson_correlation[0]])

                                        
    #print(list(zip(grid_correlation_03, correlation_value_03)))
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
