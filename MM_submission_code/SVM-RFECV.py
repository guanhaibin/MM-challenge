import csv
import numpy as np
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.datasets import make_classification

from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression

from sklearn import metrics

# Load the files as a numpy arrays. Both contain headers.
# 'Training samples' include geneIDs in col 0.
with open('DEgenes_subarray.csv', 'rb') as f:
	  reader = csv.reader(f, delimiter=',')
	  inlist = list(reader)
data = np.array(inlist)[:,:]

# 'Class labels' corresponding to 'training samples'.
with open('target_balanced.csv', 'rb') as f:
	  reader = csv.reader(f, delimiter=',')
	  inlist = list(reader)
target = np.array(inlist)[:,:]

# Remove header containing patient IDs (row 0) from both files and
# remove gene IDs (col 0) from training data file leaving only gene
# expressions.
X = data[1:,1:]
X = X.astype(np.float64)  # Convert to float64 array

y = target[1,:]
y = y.astype(np.int)  # Convert to int array

# sklearn expects X to be a matrix in which each row is one 'sample'
# and each column is one 'feature'. Therefore, need to take transpose.
X = X.T

numrows_X = len(X)
numcols_X = len(X[0])
size_y = len(y)

print "No. rows of X after transposing: %i" % numrows_X
print "No. cols of X after transposing: %i" % numcols_X
print "First 4 rows and cols of X:"
print X[0:4,0:4]
print " "
print "Size of y: %i" % size_y
print " "

# Standardize data (0 mean, 1 stdev)
scaler = StandardScaler().fit(X)
X = scaler.transform(X)

##########

# SVM-RFECV

# Support Vector Classifier (SVC) estimator.
svc = SVC(kernel="linear")

# Create the RFE object and compute a cross-validated score. 
# 'step=1' indicates that 1 feature is removed at each iteration.
# Use stratified 2-fold cross-validation for binary y. The folds
# are made by preserving the percentage of samples for each class.
# The "accuracy" scoring is proportional to the number of correct
# classifications.
rfecv = RFECV(estimator=svc, step=1, cv=StratifiedKFold(2),
              scoring='accuracy')

# Fit the RFE model and 'automatically tune' the number of selected
# features.
print "Fitting model:"
print " "
rfecv.fit(X, y)

print("Optimal number of features : %d" % rfecv.n_features_)
print " "

# Feature ranking, such that ranking_[i] corresponds to the ranking
# position of the i-th feature. Selected (i.e., estimated best)
# features are assigned rank 1.
rankings = rfecv.ranking_
print "Selector rankings:"
print rankings
print " "

# Plot number of features vs cross-validation scores
plt.figure()
plt.xlabel("Number of features selected")
plt.ylabel("Cross validation score (nb of correct classifications)")
plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)
plt.show()

# Create array to contain the location of the higest ranked features.
topranked = np.where(rankings==1)[0]
print "Position of highest-ranked features (i.e., rank=1):"
print topranked

# Extract subset of input dataset.
gene_subset = data[0,:]

# Gene positions in 'topranked' do not consider header row in 'data',
# so gene 0 in 'topranked' corresponds to gene at row 1 in 'data'.
for i in np.arange(len(topranked)):
   print topranked[i]
   gene_subset = np.vstack((gene_subset,data[topranked[i]+1,:]))

with open('SVM-RFE_subarray.csv', 'wb') as f:
    writer = csv.writer(f, delimiter=',')
    for row in gene_subset:
        writer.writerow(row)


