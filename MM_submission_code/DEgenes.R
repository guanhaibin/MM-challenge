# Load packages
library(limma)

# Print out a summary of the results
#cat("The summary of x is \n",summary(x),"\n")

# Read 'training' data. Use 'read.table' to read the data into R as a data.frame.
# Include header of patient IDs and geneIDs as row names.
data = read.table("training_balanced.csv", sep=",", header=TRUE, row.names=1)
cat("\n\n")
print("data:")
print(data[1:6,1:6]) # First 4 rows and cols
print(class(data))
cat("\n\n")

# Read 'target' vector. Include header of patient IDs.
target = read.table("target_balanced.csv", sep=",", header=TRUE, row.names=NULL)
print("target:")
print(target[1:6])
print(class(target))
cat("\n\n")

target <- unlist(target, use.names=FALSE)
print("target (after 'unlist'):")
print(target[1:6])
print(class(target))
cat("\n\n")

# Construct design matrix.
groups <- factor(c(target))
design <- model.matrix(~ groups)
print("Design matrix:")
print(head(design))
cat("\n\n")

# Fit model.
fit <- lmFit(data, design)
efit <- eBayes(fit)

# Extract significantly differentially-expressed genes.
# Interested in the second coefficient of the linear model (coeff=2),
# i.e., column 2 of 'design matrix'.
# 'genelist=NULL', since no columns of probe annotations were included.
# Filter by fold change (1x) and p-value (0.1) cutoffs.
# Adjusted using Benjimini-Hochberg false discovery rate (FDR).
# Sort by B (log-odds that the gene is differentially expressed).
DEtable <- topTable(efit, coef=2, genelist=NULL, p.value=0.1, lfc=log2(1), number=5000)

print("topTable output:")
print(head(DEtable))

# Write file of differentially-expressed genes.
write.csv(DEtable, file="DEgenes_topranked.csv")

quit()





