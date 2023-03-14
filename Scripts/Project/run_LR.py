from sklearn.model_selection import train_test_split

import numpy as np

from LR_data_preprocessing import reshape_to_1D
from LR_data_preprocessing import normalization
from LR_data_preprocessing import seperate_secondary_structures

from LR_grid_correlation import separate_bins
from LR_grid_correlation import calculate_pearson_coefficient
from LR_grid_correlation import build_high_correlation_matrix

from LR import make_Linear_Regression


def load_dataset(n_proteins):
    """522 -> small; 3374 -> big; 713 -> test"""
    with open(f"peak_matrixes_10x10_{n_proteins}.npy", "rb") as infile:
        counts_peaks_matrixes = np.load(infile)

    with open(f"Q3_percentage_{n_proteins}.npy", "rb") as infile:
        Q3_percentage = np.load(infile)

    return counts_peaks_matrixes, Q3_percentage


def build_simple_models(X_train, X_test, y_train, y_test, modelsize="big"):
    
    # train data
    helix_percentage, sheet_percentage, coil_percentage = seperate_secondary_structures(y_train)
    reshaped_matrix = reshape_to_1D(X_train)
    normalized_matrix = normalization(reshaped_matrix)


    # test data
    test_helix_percentage, test_sheet_percentage, test_coil_percentage = seperate_secondary_structures(y_test)
    test_reshaped_matrix = reshape_to_1D(X_test)
    test_normalized_matrix = normalization(test_reshaped_matrix)
    

    # Linear regression
    train_structures = [helix_percentage, sheet_percentage, coil_percentage]
    test_structures = [test_helix_percentage, test_sheet_percentage, test_coil_percentage]

    results = []
    
    for sec_structure, test_sec_structure in zip(train_structures, test_structures):
        predictions, pearson, rmse, r2, mean_predictions = make_Linear_Regression(np.array(normalized_matrix), np.array(sec_structure), np.array(test_normalized_matrix), np.array(test_sec_structure),modelsize=modelsize)
        print(sum(pearson)/len(pearson))
        results.append([predictions, pearson, rmse, r2, mean_predictions])

    return results

        
def build_rebuild_models(X_train, X_test, y_train, y_test, modelsize="big"):
    
    # train data
    helix_percentage, sheet_percentage, coil_percentage = seperate_secondary_structures(y_train)
    reshaped_matrix = reshape_to_1D(X_train)
    normalized_matrix = normalization(reshaped_matrix)

    grid_correlation_03_ = set()
    correlation_value_03_ = []

    for percentage in [helix_percentage, sheet_percentage, coil_percentage]:
        bin_by_matrix_position = separate_bins(normalized_matrix)
        grid_correlation_03, correlation_value_03 = calculate_pearson_coefficient(bin_by_matrix_position, percentage)
        grid_correlation_03_ |= set(grid_correlation_03)
        correlation_value_03_.append(grid_correlation_03)

    high_correlation_matrix = build_high_correlation_matrix(list(grid_correlation_03_), normalized_matrix)

    # test data
    test_helix_percentage, test_sheet_percentage, test_coil_percentage = seperate_secondary_structures(y_test)
    test_reshaped_matrix = reshape_to_1D(X_test)
    test_normalized_matrix = normalization(test_reshaped_matrix)

    test_high_correlation_matrix = build_high_correlation_matrix(list(grid_correlation_03_), test_normalized_matrix)

    # Linear regression
    train_structures = [helix_percentage, sheet_percentage, coil_percentage]
    test_structures = [test_helix_percentage, test_sheet_percentage, test_coil_percentage]

    results = []
    
    for sec_structure, test_sec_structure in zip(train_structures, test_structures):
        predictions, pearson, rmse, r2, mean_predictions = make_Linear_Regression(np.array(high_correlation_matrix), np.array(sec_structure), np.array(test_high_correlation_matrix), np.array(test_sec_structure),modelsize=modelsize)
        print(sum(pearson)/len(pearson))
        results.append([predictions, pearson, rmse, r2, mean_predictions])

    return results
        
    
if __name__=="__main__":

    # train data
    X_test, y_test = load_dataset(488)
    
    ## big models
    X_train, y_train = load_dataset(2353) 
    # simple big model
    res3 = build_simple_models(X_train, X_test, y_train, y_test, modelsize="big")
    # rebuild big model
    res4 = build_rebuild_models(X_train, X_test, y_train, y_test, modelsize="big")

    print()
    
    ## small models
    X_train, y_train = load_dataset(401)
    # simple small model
    res1 = build_simple_models(X_train, X_test, y_train, y_test, modelsize="small")
    # rebuild simple model
    res2 = build_rebuild_models(X_train, X_test, y_train, y_test, modelsize="small")

    print()
    names = ["helix", "sheet", "coil"]
    for i, res in enumerate([res1, res2, res3, res4]):
        print("res", i+1)
        for name, sec in zip(names, res):
            print(name)
            print("Pearson: ", round(np.mean(sec[1]),2), "+-",  round(np.std(sec[1]),2))
            print("RMSE: ", round(np.mean(sec[2]),2), "+-", round(np.std(sec[2]),2))
            print("R2: ", round(np.mean(sec[3]),2), "+-", round(np.std(sec[3]),2))
            print("---------------------------")
        print("############################")
