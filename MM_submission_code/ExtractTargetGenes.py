# Inputs training array and array that contains patient info
# (both comma-delimited). Then, creates a 'target' vector by 
# extracting only the 'Patient' IDs and 'HR_FLAG's from the
# patient info array that appear in the training array. Maintain
# the same order as in the training array. Write the 'target'
# vector to a file.
#
import csv
import numpy as np
import matplotlib.pyplot as plt
import string
from sklearn import preprocessing
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.datasets import make_classification

from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression

from sklearn import metrics

# Convert TRUE to 1 and FALSE to -1.
def to_num(s):
    return 1 if s == 'TRUE' else -1

# Compare two strings ignoring punctuation and whitespace.
def compare(s1, s2):
   remove = string.punctuation + string.whitespace
   return s1.translate(None, remove) == s2.translate(None, remove)

# File containing patient info, such as 'Patient' ID in
# the second column (col 1) and 'HR_FLAG' in the last
# column (col 43). First row (row 0) is the header.
with open('globalClinTraining.csv', 'rb') as f:
	  reader = csv.reader(f, delimiter=',')
	  inlist = list(reader)
patientInfo= np.array(inlist)[:,:]

# Dataset in which the rows are the genes and the columns
# are the patients. Cell (0,0) is empty. Row 0, beginning
# with col 1 contains 'Patient' IDs'. Column 0, beginning
# with row 1, contains 'gene' IDs. The remainder of the
# matrix contains 'gene expressions'.
with open('EMTAB4032entrezIDlevel.csv', 'rb') as f:
	  reader = csv.reader(f, delimiter=',')
	  inlist = list(reader)
geneInfo = np.array(inlist)[:,:]

numrows_patientInfo = len(patientInfo)
numcols_patientInfo = len(patientInfo[0])

numrows_geneInfo = len(geneInfo)
numcols_geneInfo = len(geneInfo[0])

print "No. rows patientInfo = %i" % numrows_patientInfo
print "No. cols patientInfo = %i" % numcols_patientInfo

print "No. rows geneInfo = %i" % numrows_geneInfo
print "No. cols geneInfo = %i" % numcols_geneInfo

# Extract 'Patient' IDs and corresponding 'HR_FLAG's.
patientFLAG = np.column_stack((patientInfo[1:,1],patientInfo[1:,43]))

# Create 'target' vector by extracting only the columns
# pertaining to 'Patient' IDs and 'HR_FLAG's.
for i in np.arange(len(patientFLAG)):
   patientFLAG[i,1] = to_num(patientFLAG[i,1])

# Form new 'target' vector that contains only the samples
# which are included in 'geneInfo' and in the same order.
cntr = 0

# Extract patient IDs from geneInfo.
geneInfoPIDs = geneInfo[0,1:]

print len(geneInfoPIDs)
print patientFLAG.shape

print patientFLAG[0,:]
print patientFLAG[0,:].shape

for i in np.arange(len(geneInfoPIDs)):
   for j in np.arange(len(patientFLAG)):
      if (compare(geneInfoPIDs[i], patientFLAG[j,0])):
         if cntr == 0:
            new_patientFLAG = patientFLAG[j,:]
            cntr += 1
         else:
            new_patientFLAG = np.vstack((new_patientFLAG,patientFLAG[j,:]))

# Take transpose for easier comparison later on.
new_patientFLAG = new_patientFLAG.T
#print new_patientFLAG

# Write 'target' vector to CSV file.
with open('target.csv', 'wb') as f:
    writer = csv.writer(f, delimiter=',')
    for row in new_patientFLAG:
        writer.writerow(row)