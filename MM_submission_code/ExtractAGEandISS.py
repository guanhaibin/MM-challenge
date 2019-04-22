# Inputs training array and array that contains patient info
# (both comma-delimited). Then, extracts the 'Patient' ID,
# 'D_Age' and 'D_ISS' from the patient info array that appear
# in the training array. Maintain the same order as in the
# training array. Write the results to a file.
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
# the second column (col 1), 'D_Age' in column 2 and
# 'D_ISS' in column 8. First row (row 0) is the header.
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

# Extract 'Patient' IDs and corresponding 'D_Age' and 'D_ISS'.
patientAgeISS = np.column_stack((patientInfo[1:,1],patientInfo[1:,2],patientInfo[1:,8]))

# Form new array that contains only the samples which
# are included in 'geneInfo' and in the same order.
cntr = 0

# Extract patient IDs from geneInfo.
geneInfoPIDs = geneInfo[0,1:]

print len(geneInfoPIDs)
print patientAgeISS.shape

print patientAgeISS[0:2,:]

for i in np.arange(len(geneInfoPIDs)):
   for j in np.arange(len(patientAgeISS)):
      if (compare(geneInfoPIDs[i], patientAgeISS[j,0])):
         if cntr == 0:
            new_patientAgeISS = patientAgeISS[j,:]
            cntr += 1
         else:
            new_patientAgeISS = np.vstack((new_patientAgeISS,patientAgeISS[j,:]))

# Take transpose for easier comparison later on.
new_patientAgeISS = new_patientAgeISS.T

# Attach fake gene IDs (888888888 to Age and 999999999 to ISS).
geneIDs = np.array((" ",888888888,999999999))
new_patientAgeISS = np.column_stack((geneIDs.T,new_patientAgeISS))

# Write results to CSV file.
with open('Age_ISS.csv', 'wb') as f:
    writer = csv.writer(f, delimiter=',')
    for row in new_patientAgeISS:
        writer.writerow(row)
