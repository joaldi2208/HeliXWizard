from sklearn.linear_model import LinearRegression
from sklearn.model_selection import RepeatedKFold
from sklearn.metrics import mean_squared_error as mae
from scipy.stats import pearsonr

import random
import numpy as np

def make_Linear_Regression(normalized_matrix, Q3_percentage, test_normalized_matrix, test_Q3_percentage, modelsize="big"):
    """makes a 10 fold cv with a linear regression model"""

    # coeffs_ = np.zeros((1,100)) # only for simple models
    r2_ = []
    pearson_corr_ = []
    rmse_ = []
    predictions_ = []

    rkf = RepeatedKFold(n_splits=10, n_repeats=100)
    for i, (train_index, test_index) in enumerate(rkf.split(Q3_percentage)):

        if modelsize == "small":
            train_index = random.sample(list(train_index), 72)

        
        reg = LinearRegression(fit_intercept=False).fit(normalized_matrix[train_index], Q3_percentage[train_index])
        # coeffs_ += reg.coef_ # only for simple models
        
        r2 = reg.score(normalized_matrix[train_index], Q3_percentage[train_index])
        predictions = reg.predict(test_normalized_matrix)

        predictions = predictions.ravel()
        predictions_.append(predictions)

        pearson_corr = pearsonr(test_Q3_percentage.ravel(), predictions)[0]
        
        pearson_corr_.append(pearson_corr)
        rmse = mae(test_Q3_percentage, predictions.ravel(), squared=False)
        rmse_.append(rmse)
        r2_.append(r2)


    mean_predictions = [sum(predictions)/len(predictions) for predictions in zip(*predictions_)]
    
    # print(coeffs_/1000) # only for simple models
    return np.concatenate(predictions_).ravel(), pearson_corr_, rmse_, r2_, mean_predictions
