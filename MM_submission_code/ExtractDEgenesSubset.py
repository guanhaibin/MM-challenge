# Minimum-redundancy maximum-revelancy (MRMR) filter.
import numpy as np
import csv
import string
from sklearn import preprocessing

# Compare two strings ignoring punctuation and whitespace.
def compare(s1, s2):
   remove = string.punctuation + string.whitespace
   return s1.translate(None, remove) == s2.translate(None, remove)

# Load the CSV files as a numpy matrices.
# Use csv reader to parse data (binary mode recommended).
# In the training file, columns are 'samples' and rows are
# 'features'.
# Read training 'data'.
with open('training_balanced.csv', 'rb') as f:
	  reader = csv.reader(f, delimiter=',')
	  inlist = list(reader)
# convert to a numpy array.
data = np.array(inlist)[:,:]

# Load CSV file containing the positions of the most
# differentially-expressed genes (rows) to be extracted.
with open('DEgenes_topranked.csv', 'rb') as f:
	  reader = csv.reader(f, delimiter=',')
	  inlist = list(reader)
# Convert to a numpy array.
topranked = np.array(inlist)

# Extract col 0 beginning with row 1
topranked = topranked[1:,0]

# Subarray of 'training_balanced.csv' containing only the most
# differentially-expressed genes. Include the header (row 0).
# Then, Append only rows pertaining to the top-ranked genes.
subarray = data[0,:]

for i in np.arange(len(topranked)):
   for j in np.arange(len(data)):
      if (compare(topranked[i], data[j,0])):
         print "i=%i  j=%i" % (i,j)
         subarray = np.vstack((subarray,data[j,:]))

with open('DEgenes_subarray.csv', 'wb') as f:
    writer = csv.writer(f, delimiter=',')
    for row in subarray:
        writer.writerow(row)
