import csv
import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import make_classification
from sklearn.decomposition import PCA
import string

from imblearn.combine import SMOTETomek

# Load the CSV files as a numpy matrices.
# Use csv reader to parse data (binary mode recommended).
# In the training file, columns are 'samples' and rows are
# 'features'.
# Read training 'data'.
with open('SVM-RFE_subarray.csv', 'rb') as f:
	  reader = csv.reader(f, delimiter=',')
	  inlist = list(reader)
# convert to a numpy array.
data = np.array(inlist)[:,:]

# Read 'Age' and 'ISS' info.
with open('AgeISS_balanced.csv', 'rb') as f:
	  reader = csv.reader(f, delimiter=',')
	  inlist = list(reader)
# convert to a numpy array.
AgeISS = np.array(inlist)[:,:]

# Remove header containing patient IDs (row 0) from training 'data'.
X = data[1:,:]

# Append 'Age' and 'ISS' features to training set.
training_AgeISS_balanced = np.vstack((AgeISS,X))

# Write the new balanced training data, file containing 'Age' and 'ISS'
# features and corresponding new 'target' vector to three separate files.
with open('training_AgeISS_balanced.csv', 'wb') as f:
    writer = csv.writer(f, delimiter=',')
    for row in training_AgeISS_balanced:
        writer.writerow(row)
