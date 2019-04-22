install.packages("missForest")
library(missForest)
##Generate 10% missing values at Random
missdta<-prodNA(test_data,noNA=0.1)

install.packages("mice")
library(mice)
imputed_Data<-mice(missdata,m=5,maxit=50,method="pmm",seed=500)

##Amelia package
install.packages("Amelia")
library(Amelia)
imputed_Data<-amelia(missdata,m=5,parallel="multicore")