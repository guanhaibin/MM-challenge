#install.packages("gbm",repos = "http://cran.us.r-project.org")
#install.packages("MLmetrics",repos = "http://cran.us.r-project.org")
#install.packages("ModelMetrics",repos = "http://cran.us.r-project.org")
#install.packages("e1071",repos = "http://cran.us.r-project.org")
#install.packages("randomForest",repos = "http://cran.us.r-project.org")
#install.packages("plyr", repos = "http://cran.us.r-project.org")
#install.packages("glmnet", repos = "http://cran.us.r-project.org")
library(glmnet)
library(gbm)
library(MLmetrics)
library(ModelMetrics)
library(caret)
library(plyr)
library(FCNN4R)
library(e1071)
library(nnet)
library(randomForest)
library(doParallel)
cl<-makePSOCKcluster(6)
registerDoParallel(cl)

##read data in R
setwd("/hguan003/MM")
DNA=read.csv("MMRF_clinical_HR_LR.attributes.33pc.csv",header=T,sep=",",row.names=1)
DNA<-as.data.frame(t(DNA))
DNA$label<-DNA[,1]
train_data<-DNA[-c(1)]
rownames(train_data)<-substr(rownames(train_data),1,9)
train_data$label<-as.factor(train_data$label)

##split into train and test data set
set.seed(123)
train.rows<-createDataPartition(y=train_data$label,p=0.7,list=FALSE)
train_MM<-train_data[train.rows,]
test_MM<-train_data[-train.rows,]
nfea<-ncol(train_MM)-1



##SVM
##the random number seed is set before each algorithm is trained to ensure that each algorithm gets the same data partitions and repeats
tc<-tune.control(cross=10)
svmtrain<-train_MM
svmfit<-svm(as.factor(label)~.,data=svmtrain,gamma=0.001,tune.control=tc)
svm_pred_train<-predict(svmfit,svmtrain[,1:nfea])
pred<-predict(svmfit,test_MM[,1:nfea])
pred<-as.vector(pred)
test_MM$label<-as.vector(test_MM$label)
SVM_accuracy<-Accuracy(pred,test_MM$label)
write.csv(pred,"SVM_pred.csv")
write.csv(SVM_accuracy,"SVM_accuracy.csv")

##random forest
set.seed(7)
fitcontrol<-trainControl(method="cv",number=10,classProbs = TRUE)
RFtrain<-train_MM
levels(RFtrain$label)<-make.names(levels(factor(RFtrain$label)))
rffit<-caret::train(label~.,
                    RFtrain,
                    method="parRF",
                    tuneGrid=expand.grid(
                      .mtry=c(30)),
                    metric="Accuracy",
                    trControl=fitcontrol)
prob_rf<-predict(rffit,test_MM[,1:nfea],type='prob')
th<-0.5
pred_rf<-factor(ifelse(prob_rf$X1>th,"1","-1"))
rf_prob_train<-predict(rffit,RFtrain[,1:nfea],type="prob")
rf_pred_train<-factor(ifelse(rf_prob_train$X1>th,"1","-1"))
prob_rf<-predict(rffit,test_MM,type='prob')
pred_rf<-as.vector(pred_rf)
test_MM$label<-as.vector(test_MM$label)
rf_accuracy<-Accuracy(pred_rf,test_MM$label)
write.csv(pred_rf,"rf_pred.csv")
write.csv(rf_accuracy,"rf_accuracy.csv")
#Learning Vector Quantization
set.seed(7)
lvqcontrol<-trainControl(method = "cv", number = 10)
lvq_grid <- expand.grid(size=c(5,10,20,50), k=c(2,4,5,7,10))
LVQtrain<-train_MM
lvq_fit<-caret::train(label~.,
                      LVQtrain,
                      method="lvq",
                      trControl=lvqcontrol,
                      tuneGrid=lvq_grid)
pred_lvq<-predict(lvq_fit,test_MM[,1:nfea])
lvq_pred_train<-predict(lvq_fit,RFtrain[,1:nfea])
pred_lvq<-as.vector(pred_lvq)
test_MM$label<-as.vector(test_MM$label)
lvq_accuracy<-Accuracy(pred_lvq,test_MM$label)
write.csv(pred_lvq,"lvq_pred.csv")
write.csv(lvq_accuracy,"lvq_accuracy.csv")
##GBM
set.seed(7)
gbm_fit<-caret::train(label~.,
                      RFtrain,
                      method="gbm",
                      verbose=FALSE,
                      trControl=fitcontrol)
prob_gbm<-predict(gbm_fit,test_MM[,1:nfea],type='prob')
th<-0.5
pred_gbm<-factor(ifelse(prob_gbm$X1>th,"1","-1"))
gbm_prob_train<-predict(gbm_fit,RFtrain[,1:nfea],type="prob")
gbm_pred_train<-factor(ifelse(gbm_prob_train$X1>th,"1","-1"))
pred_gbm<-as.vector(pred_gbm)
test_MM$label<-as.vector(test_MM$label)
gbm_accuracy<-Accuracy(pred_gbm,test_MM$label)
write.csv(gbm_accuracy,"gbm_accuracy.csv")
write.csv(gbm_pred,"gbm_pred.csv")
##glmnet Generalized Linear Model
set.seed(7)
gridsearch_for_lambda =  data.frame (alpha = 0,
                                     lambda = c (2^c(-15:15), 3^c(-15:15)))
train_control = trainControl (method="cv", number =10,
                              savePredictions =TRUE , allowParallel = FALSE,classProbs = TRUE )
glmfit<-caret::train(label~.,
                     RFtrain,
                     method="glmnet",
                     tuneGrid=gridsearch_for_lambda,
                     trControl=train_control,
                     preProcess=NULL)
prob_glm<-predict(glmfit,test_MM[,1:nfea],type='prob')
th<-0.5
pred_glm<-factor(ifelse(prob_glm$X1>th,"1","-1"))
glm_prob_train<-predict(glmfit,RFtrain[,1:nfea],type="prob")
glm_pred_train<-factor(ifelse(glm_prob_train$X1>th,"1","-1"))
pred_glm<-as.vector(pred_glm)
test_MM$label<-as.vector(test_MM$label)
glm_accuracy<-Accuracy(pred_glm,test_MM$label)
write.csv(glm_accuracy,"glm_accuracy.csv")
write.csv(pred_glm,"glm_pred.csv")


##Stacked model
#combine all the predictions of train data from above classifiers to train stacked model
train_label<-LVQtrain$label
combo_train<-data.frame(svm_pred_train,rf_pred_train,lvq_pred_train,gbm_pred_train,glm_pred_train,train_label)

for (i in 1:(ncol(combo_train)-1)){
  levels(combo_train[,i])<-levels(combo_train$train_label)
  }
fit_stacked<-caret::train(as.factor(train_label)~.,
                          combo_train,
                          method="ranger",
                          tuneGrid=expand.grid(
                            .mtry=c(2,1)),
                          metric="Accuracy",
                          trControl=lvqcontrol)
pred_stacked<-predict(fit_stacked,combo_test[,1:5])
pred_stacked<-as.vector(pred_stacked)
test_MM$label<-as.vector(test_MM$label)
stacked_accuracy<-Accuracy(pred_stacked,test_MM$label)

##output the acuuracy of all models above
ae_predict<-c(SVM_accuracy,rf_accuracy,lvq_accuracy,gbm_accuracy,glm_accuracy,stacked_accuracy)
write.csv(ae_predict,file="MMRF_accuracy.csv")

#save model seperately
saveRDS(svmfit,"./M_svmfit.rds")
saveRDS(rffit,"./M_rffit.rds")
saveRDS(lvq_fit,"./M_lvq_fit.rds")
saveRDS(gbm_fit,"./M_gbm_fit.rds")
saveRDS(glmfit,"./M_glmfit.rds")
saveRDS(fit_stacked,"./M_fit_stacked.rds")


stopImplicitCluster()
