library(RSQLite)
library(plyr)
library(FCNN4R)
library(MLmetrics)
library(ModelMetrics)
library(stringr)
library(e1071)
library(caret)
library(doParallel)
cl<-makePSOCKcluster(6)
registerDoParallel(cl)
rm(list=ls())
setwd("U://hguan003/MM")
U_genes_list<-read.csv("U_genes.csv",row.names = 1)
U_genes_list<-as.character(U_genes_list[,1])
E_genes_list<-read.csv("E_genes.csv",row.names = 1)
E_genes_list<-as.character(E_genes_list[,1])
H_genes_list<-read.csv("H_genes.csv",row.names = 1)
H_genes_list<-as.character(H_genes_list[,1])

Gene_expression=read.csv("HOVON_subarray.csv",header=T,sep=",",row.names=1)
Gene_expression<-as.data.frame(t(scale(Gene_expression)))
HR_FLAG<-read.csv("HOVON_target.csv",header=T,sep=",")
HR_FLAG<-t(HR_FLAG)
train_data<-merge(Gene_expression,HR_FLAG,by="row.names",all.x=TRUE)
row.names(train_data)<-train_data[,1]
train_data<-train_data[-c(1)]
colnames(train_data)[ncol(train_data)]<-c("label")
train_data$label<-as.factor(train_data$label)
MMRF_test<-train_data

#split MMRF test data into different subset to fit in different models
MMRF_U_test<-setNames(as.data.frame(matrix(0,nrow=nrow(MMRF_test),ncol=length(U_genes_list))),U_genes_list)
idx<-colnames(MMRF_test)[colnames(MMRF_test)%in%colnames(MMRF_U_test)]
MMRF_U_test[,idx]<-MMRF_test[,idx]
rownames(MMRF_U_test)<-rownames(MMRF_test)
MMRF_U_test$label<-MMRF_test$label
MMRF_E_test<-setNames(as.data.frame(matrix(0,nrow=nrow(MMRF_test),ncol=length(E_genes_list))),E_genes_list)
idx<-colnames(MMRF_test)[colnames(MMRF_test)%in%colnames(MMRF_E_test)]
MMRF_E_test[,idx]<-MMRF_test[,idx]
rownames(MMRF_E_test)<-rownames(MMRF_test)
MMRF_E_test$label<-MMRF_test$label
MMRF_H_test<-setNames(as.data.frame(matrix(0,nrow=nrow(MMRF_test),ncol=length(H_genes_list))),H_genes_list)
idx<-colnames(MMRF_test)[colnames(MMRF_test)%in%colnames(MMRF_H_test)]
MMRF_H_test[,idx]<-MMRF_test[,idx]
rownames(MMRF_H_test)<-rownames(MMRF_test)
MMRF_H_test$label<-MMRF_test$label



##combine all test data correspomding to the different model
test_U<-MMRF_U_test
test_E<-MMRF_E_test
test_H<-MMRF_H_test

th<-0.5

##load the U_mdoel 
U_svmfit<-readRDS("./U_svmfit.rds")
test_U_svm_prob<-predict(U_svmfit,test_U[,1:(ncol(test_U)-1)],probability = TRUE)
test_U_svm<-predict(U_svmfit,test_U[,1:(ncol(test_U)-1)])
U_rffit<-readRDS("./U_rffit.rds")
test_U_rf_prob<-predict(U_rffit,test_U[,1:(ncol(test_U)-1)],type="prob")
test_U_rf<-factor(ifelse(test_U_rf_prob$X1>th,"1","-1"))
U_nn_fit<-readRDS("./U_nn_fit.rds")
test_U_nn_prob<-predict(U_nn_fit,test_U[,1:(ncol(test_U)-1)],type="prob")
test_U_nn<-predict(U_nn_fit,test_U[,1:(ncol(test_U)-1)])
U_lvq_fit<-readRDS("./U_lvq_fit.rds")
test_U_lvq<-predict(U_lvq_fit,test_U[,1:(ncol(test_U)-1)])
U_gbm_fit<-readRDS("./U_gbm_fit.rds")
test_U_gbm_prob<-predict(U_gbm_fit,test_U[,1:(ncol(test_U)-1)],type="prob")
test_U_gbm<-factor(ifelse(test_U_gbm_prob$X1>th,"1","-1"))
U_glmfit<-readRDS("./U_glmfit.rds")
test_U_glm_prob<-predict(U_glmfit,test_U[,1:(ncol(test_U)-1)],type="prob")
test_U_glm<-factor(ifelse(test_U_glm_prob$X1>th,"1","-1"))
U_fit_stacked<-readRDS("./U_fit_stacked.rds")
combo_U_test<-data.frame(test_U_svm,test_U_rf,test_U_nn,test_U_lvq,test_U_gbm,test_U_glm,label<-as.factor(test_U$label))
colnames(combo_U_test)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train","train_label" )
test_U_prob<-predict(U_fit_stacked,combo_U_test,type="prob")
th<-0.5
test_U_stacked<-factor(ifelse(test_U_prob$X1>th,"1","-1"))


##load the E_mdoel
E_svmfit<-readRDS("./E_svmfit.rds")
test_E_svm_prob<-predict(E_svmfit,test_E[,1:(ncol(test_E)-1)],probability = TRUE)
test_E_svm<-predict(E_svmfit,test_E[,1:(ncol(test_E)-1)])
E_rffit<-readRDS("./E_rffit.rds")
test_E_rf_prob<-predict(E_rffit,test_E[,1:(ncol(test_E)-1)],type="prob")
test_E_rf<-factor(ifelse(test_E_rf_prob$X1>th,"1","-1"))
E_nn_fit<-readRDS("./E_nn_fit.rds")
test_E_nn_prob<-predict(E_nn_fit,test_E[,1:(ncol(test_E)-1)],type="prob")
test_E_nn<-predict(E_nn_fit,test_E[,1:(ncol(test_E)-1)])
E_lvq_fit<-readRDS("./E_lvq_fit.rds")
test_E_lvq<-predict(E_lvq_fit,test_E[,1:(ncol(test_E)-1)])
E_gbm_fit<-readRDS("./E_gbm_fit.rds")
test_E_gbm_prob<-predict(E_gbm_fit,test_E[,1:(ncol(test_E)-1)],type="prob")
test_E_gbm<-factor(ifelse(test_E_gbm_prob$X1>th,"1","-1"))
E_glmfit<-readRDS("./E_glmfit.rds")
test_E_glm_prob<-predict(E_glmfit,test_E[,1:(ncol(test_E)-1)],type="prob")
test_E_glm<-factor(ifelse(test_E_glm_prob$X1>th,"1","-1"))
E_fit_stacked<-readRDS("./E_fit_stacked.rds")
combo_E_test<-data.frame(test_E_svm,test_E_rf,test_E_nn,test_E_lvq,test_E_gbm,test_E_glm,label<-as.factor(test_E$label))
colnames(combo_E_test)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train" ,"train_label" )
test_E_prob<-predict(E_fit_stacked,combo_E_test,type="prob")
th<-0.5
test_E_stacked<-factor(ifelse(test_E_prob$X1>th,"1","-1"))

##load the H_mdoel
H_svmfit<-readRDS("./H_svmfit.rds")
test_H_svm_prob<-predict(H_svmfit,test_H[,1:(ncol(test_H)-1)],probability = TRUE)
test_H_svm<-predict(H_svmfit,test_H[,1:(ncol(test_H)-1)])
H_rffit<-readRDS("./H_rffit.rds")
test_H_rf_prob<-predict(H_rffit,test_H[,1:(ncol(test_H)-1)],type="prob")
test_H_rf<-factor(ifelse(test_H_rf_prob$X1>th,"1","-1"))
H_nn_fit<-readRDS("./H_nn_fit.rds")
test_H_nn_prob<-predict(H_nn_fit,test_H[,1:(ncol(test_H)-1)],type="prob")
test_H_nn<-predict(H_nn_fit,test_H[,1:(ncol(test_H)-1)])
H_lvq_fit<-readRDS("./H_lvq_fit.rds")
test_H_lvq<-predict(H_lvq_fit,test_H[,1:(ncol(test_H)-1)])
H_gbm_fit<-readRDS("./H_gbm_fit.rds")
test_H_gbm_prob<-predict(H_gbm_fit,test_H[,1:(ncol(test_H)-1)],type="prob")
test_H_gbm<-factor(ifelse(test_H_gbm_prob$X1>th,"1","-1"))
H_glmfit<-readRDS("./H_glmfit.rds")
test_H_glm_prob<-predict(H_glmfit,test_H[,1:(ncol(test_H)-1)],type="prob")
test_H_glm<-factor(ifelse(test_H_glm_prob$X1>th,"1","-1"))
H_fit_stacked<-readRDS("./H_fit_stacked.rds")
combo_H_test<-data.frame(test_H_svm,test_H_rf,test_H_nn,test_H_lvq,test_H_gbm,test_H_glm,label<-as.factor(test_H$label))
colnames(combo_H_test)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train","train_label" )
test_H_prob<-predict(H_fit_stacked,combo_H_test,type="prob")
th<-0.5
test_H_stacked<-factor(ifelse(test_H_prob$X1>th,"1","-1"))

##load ensemble model
ensemble_rf<-readRDS("./ensemble_rf.rds")
ensemble_nn<-readRDS("./ensemble_nn.rds")
##global ensemble_test data
U_models_test<-cbind(combo_U_test[,1:6],test_U_stacked)
colnames(U_models_test)<-gsub("_pred_","_U_",colnames(U_models_test))
E_models_test<-cbind(combo_E_test[,1:6],test_E_stacked)
colnames(E_models_test)<-gsub("_pred_","_E_",colnames(E_models_test))
H_models_test<-cbind(combo_H_test[,1:6],test_H_stacked)
colnames(H_models_test)<-gsub("_pred_","_H_",colnames(H_models_test))
ensemble_test<-cbind(U_models_test,E_models_test,H_models_test)
ensemble_test$label<-combo_U_test$train_label
colnames(ensemble_test)[1:21]<-gsub("1","",ensemble_rf$coefnames)


ensemble_rf_prob<-predict(ensemble_rf,ensemble_test[,1:21],type="prob")  ###
th<-0.5
ensemble_rf_pred<-factor(ifelse(ensemble_rf_prob$X1>th,"1","-1"))
ensemble_rf_accuracy<-Accuracy(ensemble_rf_pred,ensemble_test$label)
ensemble_nn_prob<-predict(ensemble_nn,ensemble_test[,1:21],type="prob")
th<-0.5
ensemble_nn_pred<-factor(ifelse(ensemble_nn_prob$X1>th,"1","-1"))
ensemble_nn_accuracy<-Accuracy(ensemble_nn_pred,ensemble_test$label)

ensemble_accuracy<-as.data.frame(matrix(0,ncol=1,nrow = 23))
for (i in 1:(ncol(ensemble_test)-1)){
  levels(ensemble_test[,i])<-levels(ensemble_test$label)
  a<-Accuracy(ensemble_test[,i],ensemble_test$label)
  ensemble_accuracy[i,1]<-a
  rownames(ensemble_accuracy)[i]<-colnames(ensemble_test)[i]
}

ensemble_accuracy[22,1]<-ensemble_rf_accuracy
rownames(ensemble_accuracy)[22]<-"ensemble_rf"
ensemble_accuracy[23,1]<-ensemble_nn_accuracy
rownames(ensemble_accuracy)[23]<-"ensemble_nn"
rownames(ensemble_accuracy)<-gsub("train","test",rownames(ensemble_accuracy))
write.csv(ensemble_accuracy,"MMRF_ensemble_accuracy.csv")



###generate prediction score
predictionscore<-array()
PS_U<-array()
PS_H<-array()
PS_E<-array()
for (i in 1:length(ensemble_rf_pred)){
  PS_U[i]<-sum(attr(test_U_svm_prob,"probabilities")[i,2],test_U_gbm_prob$X1[i],test_U_glm_prob$X1[i],test_U_nn_prob[i,"1"],test_U_rf_prob[i,"1"],test_U_prob$X1[i])
  PS_H[i]<-sum(attr(test_H_svm_prob,"probabilities")[i,2],test_H_gbm_prob$X1[i],test_H_glm_prob$X1[i],test_H_nn_prob[i,"1"],test_H_rf_prob[i,"1"],test_H_prob$X1[i])
  PS_E[i]<-sum(attr(test_E_svm_prob,"probabilities")[i,2],test_E_gbm_prob$X1[i],test_E_glm_prob$X1[i],test_E_nn_prob[i,"1"],test_E_rf_prob[i,"1"],test_E_prob$X1[i])
  if (ensemble_rf_pred[i]=="1"){
    predictionscore[i]<-max(c(PS_U[i],PS_H[i],PS_E[i]))
  } else {
    predictionscore[i]<-min(c(PS_U[i],PS_H[i],PS_E[i]))
  }
}




predictionscore<-round(predictionscore,5)
highriskflag<-ensemble_rf_pred
Patient<-rownames(ensemble_test)


###Output1: max P(E1),P(E2),P(E3)
output<-data.frame(Patient,predictionscore,highriskflag)
write.csv(output,"MMRF_predictionfile.csv",row.names = FALSE)



stopImplicitCluster()



