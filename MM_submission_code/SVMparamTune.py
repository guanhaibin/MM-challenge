# The SVM classifier is optimized by cross-validation using
# GridSearchCV. The data is split into a 'development' set
# (to be fed to the GridSearchCV instance) for training the
# 'C' and 'gamma' hyper-parameters of an 'rbf-kernel' SVM
# and an 'evaluation' set to compute performance metrics.

import numpy as np
import csv
import pylab as pl

from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score
from sklearn import svm
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.metrics import classification_report

# Load the files as a numpy arrays. Both contain headers.
# 'Training samples' include geneIDs in col 0.
with open('training_AgeISS_balanced.csv', 'rb') as f:
	  reader = csv.reader(f, delimiter=',')
	  inlist = list(reader)
data = np.array(inlist)[:,:]

# 'Class labels' corresponding to 'training samples'.
with open('target_balanced.csv', 'rb') as f:
	  reader = csv.reader(f, delimiter=',')
	  inlist = list(reader)
target = np.array(inlist)[:,:]

# 'sklearn' expects X to be a matrix in which every row is one
# 'sample' and every column is one 'feature'. Therefore, need
# to take transpose.
X = data[1:,1:] # remove header and geneIDs
X = X.astype(np.float64)
X = X.T
print "X:" 
print X.shape
print(" ")

y = target[1:,:] # remove header
y = y.astype(np.int)  # convert to 'int'
y = np.ravel(y) # convert to 1-dim array
print "y:"
print len(y)
print(" ")

# It is usually a good idea to scale the data for SVM training.
# We are cheating a bit in this example in scaling all of the data,
# instead of fitting the transformation on the trainingset and
# just applying it on the test set.
# Standardize data (0 mean, 1 stdev)
scaler = StandardScaler().fit(X)
X = scaler.transform(X)

##########

# Split data so that 60% is used for training while leaving 40% for
# testing (evaluating) the classifier.
# By default, the data is 'shuffled' before the split.
# The 'stratify' parameter is used to maintain the same proportion of each
# class when splitting (i.e., it will make sure that each subset resulting
# from the split has the same percentage of 1 and -1 y values as the full
# dataset).
# Set 'random_state' to an integer to get identical results for each split.
# Otherwise, the shuffling will be different each time 'KFold' is iterated.
#X_train, X_test, y_train, y_test = train_test_split(X, y,
#                          test_size=0.4, random_state=0)
X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y,
                                       test_size=0.4, random_state=0)

# Set the parameters by cross-validation
tuned_parameters = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4],
                     'C': [1, 10, 100, 1000]},
                    {'kernel': ['linear'], 'C': [1, 10, 100, 1000]}]

scores = ['precision', 'recall']

for score in scores:
    print("Tuning hyper-parameters for %s" % score)
    print(" ")

    # Use 5-fold cross validation to tune parameters.
    clf = GridSearchCV(SVC(), tuned_parameters, cv=5,
                       scoring='%s_macro' % score)
    clf.fit(X_train, y_train)

    print("Best parameters set found on development set:")
    print(" ")
    print(clf.best_params_)
    print(" ")
    print("Grid scores on development set:")
    print(" ")
    means = clf.cv_results_['mean_test_score']
    stds = clf.cv_results_['std_test_score']
    for mean, std, params in zip(means, stds, clf.cv_results_['params']):
        print("%0.3f (+/-%0.03f) for %r"
              % (mean, std * 2, params))
    print(" ")

    print("Detailed classification report:")
    print(" ")
    print("The model is trained on the full development set.")
    print("The scores are computed on the full evaluation set.")
    print(" ")
    y_true, y_pred = y_test, clf.predict(X_test)
    print(classification_report(y_true, y_pred))
    print(" ")

