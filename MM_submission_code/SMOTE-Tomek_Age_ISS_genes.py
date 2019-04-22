import csv
import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import make_classification
from sklearn.decomposition import PCA
import string

from imblearn.combine import SMOTETomek

# Compare two strings ignoring punctuation and whitespace.
def compare(s1, s2):
   remove = string.punctuation + string.whitespace
   return s1.translate(None, remove) == s2.translate(None, remove)

# Load the CSV files as a numpy matrices.
# Use csv reader to parse data (binary mode recommended).
# In the training file, columns are 'samples' and rows are
# 'features'.
# Read training 'data'.
with open('EMTAB4032entrezIDlevel.csv', 'rb') as f:
	  reader = csv.reader(f, delimiter=',')
	  inlist = list(reader)
# convert to a numpy array.
data = np.array(inlist)[:,:]

# Read 'target' vector.
with open('target.csv', 'rb') as f:
	  reader = csv.reader(f, delimiter=',')
	  inlist = list(reader)
# convert to a numpy array.
target = np.array(inlist)[:,:]

# Read 'Age' and 'ISS' info.
with open('Age_ISS.csv', 'rb') as f:
	  reader = csv.reader(f, delimiter=',')
	  inlist = list(reader)
# convert to a numpy array.
AgeISS = np.array(inlist)[:,:]

# Remove header containing patient IDs (row 0) from all files and
# remove gene IDs (col 0) from training data and AgeISS files.
y = target[1,:]
size_y = len(y)
print "length y: %i" % size_y

X = data[1:,1:]
numcols_X = len(X[0])
numrows_X = len(X)
print "numcols X: %i" % numcols_X
print "numrows X: %i" % numrows_X

AgeISSinfo = AgeISS[1:,1:]
numcols_AgeISSinfo = len(AgeISSinfo[0])
numrows_AgeISSinfo = len(AgeISSinfo)
print "numcols AgeISSinfo: %i" % numcols_AgeISSinfo
print "numrows AgeISSinfo: %i" % numrows_AgeISSinfo
print " "

# Attach 'Age' and 'ISS' features to training set.
X = np.vstack((AgeISSinfo,X))

print "X with Age and ISS:"
print X[:3,:3]
print " "

cntr = 0

# Remove samples from the dataset that either have an 'Age' or an 'ISS' value of NA.
for j in np.arange(len(X[0])):
   if (not (compare("NA", X[1,j]) or compare("NA", X[2,j]))):
      if cntr == 0:
         new_X = X[:,j]
         new_y = y[j]
         cntr += 1
      else:
         new_X = np.column_stack((new_X,X[:,j]))
         new_y = np.append(new_y,y[j])

X = new_X
X = X.astype(np.float64)  # Convert to float64 array

y = new_y
y = y.astype(np.int)  # Convert to int array

size_y = len(y)
print "length y: %i" % size_y
numcols_X = len(X[0])
numrows_X = len(X)
print "numcols X: %i" % numcols_X
print "numrows X: %i" % numrows_X
print " "

# sklearn expects X to be a matrix in which each row is one 'sample'
# and each column is one 'feature'. Therefore, need to take transpose.
X = X.T

numcols_X = len(X[0])
numrows_X = len(X)
size_y = len(y)

print "cols X (after transposing): %i" % numcols_X
print "rows X (after transposing): %i" % numrows_X
print " "
print "size y: %i" % size_y
print "y:"
print y
print " "

cntneg1=0
cntpos1=0

for i in np.arange(size_y):
    if y[i] < 0:
        cntneg1+=1
    else:
        cntpos1+=1

print "No. of -1's: %i" % cntneg1
print "No. of 1's: %i" % cntpos1
print " "

# Apply SMOTE + Tomek links
sm = SMOTETomek()
X_balanced, y_balanced = sm.fit_sample(X, y)

numcols_X_balanced = len(X_balanced[0])
numrows_X_balanced = len(X_balanced)
size_y_balanced = len(y_balanced)

print "cols X_balanced: %i" % numcols_X_balanced
print "rows X_balanced: %i" % numrows_X_balanced
print " "
print "size y_balanced: %i" % size_y_balanced
print "y_balanced:"
print y_balanced
print " "

cntneg1=0
cntpos1=0

for i in np.arange(size_y_balanced):
    if y_balanced[i] < 0:
        cntneg1+=1
    else:
        cntpos1+=1

print "No. of -1's: %i" % cntneg1
print "No. of 1's: %i" % cntpos1
print " "

print "Increase in no. of new samples after balancing: %i" % (size_y_balanced - size_y)
print " "

# The 'target' vector contains new patient IDs for the balanced dataset
# (artificially generated) in row 0 and corresponding new HR_FLAGs in row 1.
target_balanced_IDs = np.arange(size_y_balanced)
target_balanced = np.vstack((target_balanced_IDs,y_balanced))

# Take transpose of 'X_balanced' to match original dataset.
X_balanced = X_balanced.T

print "X_balanced:"
print X_balanced[:5,:3]
print " "

# Separate 'Age' and 'ISS' features.
AgeISSinfo_balanced = X_balanced[:2,:]
X_balanced = X_balanced[2:,:]

print "AgeISSinfo_balanced:"
print AgeISSinfo_balanced[:,:3]
print " "

print "X_balanced (after separationg 'Age' and 'ISS'):"
print X_balanced[:3,:3]
print " "

# Add gene IDs to 'AgeISSinfo' and 'X_balanced'.
data_balanced = np.column_stack((data[1:,0],X_balanced))
print "data_balanced:"
print data_balanced[:3,:3]
print " "

AgeISSinfo_balanced = np.column_stack((AgeISS[1:,0],AgeISSinfo_balanced))
print "AgeISSinfo_balanced:"
print AgeISSinfo_balanced[:3,:3]
print " "

# Then, paste header containing new patient IDs, beginning with a 'blank'
# to all three.
header = [" "]
print "header:"
print header
header = np.append(header,target_balanced_IDs)
print "header:"
print header
print(header.shape)
data_balanced = np.vstack((header,data_balanced))
AgeISSinfo_balanced = np.vstack((header,AgeISSinfo_balanced))

print "No. rows of data_balanced (incl. patient IDs): %i" % len(data_balanced)
print "No. cols of data_balanced (incl. gene IDs): %i" % len(data_balanced[0])
print "No. rows of AgeISSinfo_balanced (incl. patient IDs): %i" % len(AgeISSinfo_balanced)
print "No. cols of AgeISSinfo_balanced (incl. gene IDs): %i" % len(AgeISSinfo_balanced[0])

# Write the new balanced training data, file containing 'Age' and 'ISS'
# features and corresponding new 'target' vector to three separate files.
with open('training_balanced.csv', 'wb') as f:
    writer = csv.writer(f, delimiter=',')
    for row in data_balanced:
        writer.writerow(row)

with open('AgeISS_balanced.csv', 'wb') as f:
    writer = csv.writer(f, delimiter=',')
    for row in AgeISSinfo_balanced:
        writer.writerow(row)

with open('target_balanced.csv', 'wb') as f:
    writer = csv.writer(f, delimiter=',')
    for row in target_balanced:
        writer.writerow(row)
