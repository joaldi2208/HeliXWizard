import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split

from sklearn.multioutput import MultiOutputRegressor
from sklearn.linear_model import Ridge

from sklearn.ensemble import RandomForestRegressor

with open("peak_matrixes_10x10.npy", "rb") as infile:
    counts_peaks_matrixes = np.load(infile)
    print(counts_peaks_matrixes[0].shape)

with open("Q3_percentage.npy", "rb") as infile:
    Q3_percentage = np.load(infile)

X_train, X_test, y_train, y_test = train_test_split(counts_peaks_matrixes, Q3_percentage, test_size=0.25, random_state=73)


#print(X_train[0])
# --> C E H 
X_train = np.reshape(X_train, (len(X_train),100))
X_test = np.reshape(X_test, (len(X_test),100)) # das muss noch normiert werden
#print(X_train[0])

# normalization
normalized_X_train = []
for hsqc in X_train:
    total_count = np.sum(hsqc)
    normalized_hsqc = []
    for grid in hsqc:
        normalized_grid = grid / total_count
        normalized_hsqc.append(normalized_grid)
    normalized_X_train.append(normalized_hsqc)

# training data sets
y_train_helix = [[sec[2]] for sec in y_train]
y_test_helix = [[sec[2]] for sec in y_test]

y_train_sheet = [[sec[1]] for sec in y_train]
y_test_sheet = [[sec[1]] for sec in y_test]

y_train_coil = [[sec[0]] for sec in y_train]
y_test_coil = [[sec[0]] for sec in y_test]


# fitting and predicting helix
print("Helix")
reg = LinearRegression().fit(normalized_X_train, y_train_helix)
preds = reg.predict(normalized_X_train).ravel()

y_test_helix = np.array(y_test_helix).ravel()

print(pearsonr(preds, y_train_helix))
print("--------------------------------")
plt.scatter(preds, y_train_helix)
plt.show()

# fitting and predicting sheet
print("sheet")
reg = LinearRegression().fit(normalized_X_train, y_train_sheet)
preds = reg.predict(normalized_X_train).ravel()

y_test_sheet = np.array(y_test_sheet).ravel()

print(pearsonr(preds, y_train_sheet))
print("--------------------------------")
plt.scatter(y_train_sheet, preds)
plt.show()

# fitting and predicting coil
print("Coil")
reg = LinearRegression().fit(normalized_X_train, y_train_coil)
preds = reg.predict(normalized_X_train).ravel()

y_test_coil = np.array(y_test_sheet).ravel()

print(pearsonr(preds, y_train_coil))
print("--------------------------------")
plt.scatter(preds, y_train_coil)
plt.show()

#rf = RandomForestRegressor().fit(X_train, y_train)
#preds = rf.predict(X_test).ravel()

#print(preds[:2])
#print(y_test[:2])
#print(pearsonr(preds, y_test))

#import matplotlib.pyplot as plt

#plt.scatter(preds, y_test)
#plt.show()
