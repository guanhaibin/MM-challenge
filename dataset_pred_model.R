library(glmnet)
library(xgboost)
library(plyr)
library(FCNN4R)
library(MLmetrics)
library(ModelMetrics)
library(stringr)
library(e1071)
library(caret)
library(matrixStats)
library(Linnorm)


setwd("U:/hguan003/MM_Resample")
U_genes_list<-read.csv("U_genes.csv",row.names = 1)
U_genes_list<-as.character(U_genes_list[,1])
E_genes_list<-read.csv("E_genes.csv",row.names = 1)
E_genes_list<-as.character(E_genes_list[,1])
H_genes_list<-read.csv("H_genes.csv",row.names = 1)
H_genes_list<-as.character(H_genes_list[,1])
M_genes_list<-read.csv("M_genes.csv",row.names = 1)
M_genes_list<-as.character(M_genes_list[,1])

U_genes_eu_list<-read.csv("U_genes_eu.csv",row.names = 1)
U_genes_eu_list<-as.character(U_genes_eu_list[,1])
E_genes_eu_list<-read.csv("E_genes_eu.csv",row.names = 1)
E_genes_eu_list<-as.character(E_genes_eu_list[,1])
M_genes_eu_list<-read.csv("M_genes_eu.csv",row.names = 1)
M_genes_eu_list<-as.character(M_genes_eu_list[,1])


load("emc92_uams70.Rd")
genes_list<-unique(c(emc92,uams70))
eu_genes_list<-unique(genes_list)
model.entrezs<-c(U_genes_eu_list,E_genes_eu_list,M_genes_eu_list)
model.entrezs<-unique(model.entrezs)


#read UAMS datasets with selected genes
Gene_expression=read.csv("UAMS_balanced_Ensemble.csv",header=T,sep=",",row.names=1)
Gene_expression<-as.data.frame(t(Gene_expression))
HR_FLAG<-read.csv("UAMS_target.csv",header=T,sep=",")
HR_FLAG<-t(HR_FLAG)
train_data<-merge(Gene_expression,HR_FLAG,by="row.names",all.x=TRUE)
row.names(train_data)<-train_data[,1]
train_data<-train_data[-c(1)]
colnames(train_data)[ncol(train_data)]<-c("label")
train_data$label<-as.factor(train_data$label)
rownames(train_data)<-gsub("X","U",rownames(train_data))
UAMS_label<-train_data$label

select.data<-setNames(as.data.frame(matrix(0,nrow=nrow(train_data),ncol=length(model.entrezs))),model.entrezs)
idx<-colnames(train_data)[colnames(train_data)%in%colnames(select.data)]
select.data[,idx]<-train_data[,idx]
rownames(select.data)<-rownames(train_data)
train_data <- select.data
#train_data$set<-"UAMS"
train_data$set<-"1"
train_data$label<-UAMS_label

U_test_MM<-read.csv("U_test_MM.csv",row.names = 1)
U_train_MM<-read.csv("U_train_MM.csv",row.names = 1)
U_test_MM<-as.character(U_test_MM[,1])
U_train_MM<-as.character(U_train_MM[,1])
UAMS_train<-subset(train_data,rownames(train_data)%in%U_train_MM)
UAMS_test<-subset(train_data,rownames(train_data)%in%U_test_MM)


#EMTAB
Gene_expression=read.csv("EMTAB_training_balanced.csv",header=T,sep=",",row.names=1)
Gene_expression<-as.data.frame(t(Gene_expression))
HR_FLAG<-read.csv("EMTAB_target_balanced.csv",header=T,sep=",")
HR_FLAG<-t(HR_FLAG)
train_data<-merge(Gene_expression,HR_FLAG,by="row.names",all.x=TRUE)
row.names(train_data)<-train_data[,1]
train_data<-train_data[-c(1)]
colnames(train_data)[ncol(train_data)]<-c("label")
train_data$label<-as.factor(train_data$label)
rownames(train_data)<-gsub("X","E",rownames(train_data))
EMTAB_label<-train_data$label


select.data<-setNames(as.data.frame(matrix(0,nrow=nrow(train_data),ncol=length(model.entrezs))),model.entrezs)
idx<-colnames(train_data)[colnames(train_data)%in%colnames(select.data)]
select.data[,idx]<-train_data[,idx]
rownames(select.data)<-rownames(train_data)
train_data <- select.data
#train_data$set<-"EMTAB"
train_data$set<-"2"
train_data$label<-EMTAB_label


E_train_MM<-read.csv("E_train_MM.csv",row.names = 1)
E_train_MM<-as.character(E_train_MM[,1])
E_test_MM<-read.csv("E_test_MM.csv",row.names = 1)
E_test_MM<-as.character(E_test_MM[,1])
EMTAB_train<-subset(train_data,rownames(train_data)%in%E_train_MM)
EMTAB_test<-subset(train_data,rownames(train_data)%in%E_test_MM)

#read HOVON datasets with selected genes
Gene_expression=read.csv("HOVON_training_balanced.csv",header=T,sep=",",row.names=1)
Gene_expression<-as.data.frame(t(Gene_expression))
HR_FLAG<-read.csv("HOVON_target_balanced.csv",header=T,sep=",")
HR_FLAG<-t(HR_FLAG)
train_data<-merge(Gene_expression,HR_FLAG,by="row.names",all.x=TRUE)
row.names(train_data)<-train_data[,1]
train_data<-train_data[-c(1)]
colnames(train_data)[ncol(train_data)]<-c("label")
train_data$label<-as.factor(train_data$label)
rownames(train_data)<-gsub("X","H",rownames(train_data))
HOVON_label<-train_data$label

select.data<-setNames(as.data.frame(matrix(0,nrow=nrow(train_data),ncol=length(model.entrezs))),model.entrezs)
idx<-colnames(train_data)[colnames(train_data)%in%colnames(select.data)]
select.data[,idx]<-train_data[,idx]
rownames(select.data)<-rownames(train_data)
train_data <- select.data
#train_data$set<-"HOVON"
train_data$set<-"1"
train_data$label<-HOVON_label

H_train_MM<-read.csv("H_train_MM.csv",row.names = 1)
H_train_MM<-as.character(H_train_MM[,1])
H_test_MM<-read.csv("H_test_MM.csv",row.names = 1)
H_test_MM<-as.character(H_test_MM[,1])
HOVON_train<-subset(train_data,rownames(train_data)%in%H_train_MM)
HOVON_test<-subset(train_data,rownames(train_data)%in%H_test_MM)

#read MMRF datasets with selected genes
Gene_expression=read.csv("MMRF_training_balanced.csv",header=T,sep=",",row.names=1)
#Gene_expression<-as.data.frame(t(scale(Gene_expression)))
Gene_expression<-as.data.frame(t(Gene_expression))
HR_FLAG<-read.csv("MMRF_target.csv",header=T,sep=",")
HR_FLAG<-t(HR_FLAG)
train_data<-merge(Gene_expression,HR_FLAG,by="row.names",all.x=TRUE)
row.names(train_data)<-train_data[,1]
train_data<-train_data[-c(1)]
colnames(train_data)[ncol(train_data)]<-c("label")
train_data$label<-as.factor(train_data$label)
rownames(train_data)<-gsub("X","M",rownames(train_data))
MMRF_label<-train_data$label

select.data<-setNames(as.data.frame(matrix(0,nrow=nrow(train_data),ncol=length(model.entrezs))),model.entrezs)
idx<-colnames(train_data)[colnames(train_data)%in%colnames(select.data)]
select.data[,idx]<-train_data[,idx]
rownames(select.data)<-rownames(train_data)
train_data <- select.data
#train_data$set<-"MMRF"
train_data$set<-"3"
train_data$label<-MMRF_label

M_train_MM<-read.csv("M_train_MM.csv",row.names = 1)
M_train_MM<-as.character(M_train_MM[,1])
M_test_MM<-read.csv("M_test_MM.csv",row.names = 1)
M_test_MM<-as.character(M_test_MM[,1])
MMRF_train<-subset(train_data,rownames(train_data)%in%M_train_MM)
MMRF_test<-subset(train_data,rownames(train_data)%in%M_test_MM)


train_data<-as.data.frame(rbind(UAMS_train,EMTAB_train,HOVON_train,MMRF_train))
test_data<-as.data.frame(rbind(UAMS_test,EMTAB_test,HOVON_test,MMRF_test))
train_MM<-train_data
test_MM<-test_data
train_MM_label<-train_MM$label
train_MM$label<-NULL
test_MM_label<-test_MM$label
test_MM$label<-NULL
train_MM_set<-train_MM$set
train_MM$set<-NULL
test_MM_set<-test_MM$set
test_MM$set<-NULL

###scale test and train set with trai_mean and train_sd
set_mean<-rowMeans(as.matrix(t(train_MM)))
set_sd<-rowSds(as.matrix(t(train_MM)))


train=sweep(train_MM,MARGIN=2,set_mean,"-")
train_MM<-sweep(train,MARGIN=2,set_sd,"/")
test=sweep(test_MM,MARGIN=2,set_mean,"-")
test_MM<-sweep(test,MARGIN=2,set_sd,"/")
train_MM$set<-train_MM_set
test_MM$set<-test_MM_set
nfea<-ncol(train_MM)-1


##create a dataframe to store prediction score for different model
ae_predict<-setNames(as.data.frame(matrix(ncol=3,nrow = 0)),c("model","Accuracy","classloss"))
class_result<-as.data.frame(matrix(ncol=nrow(test_MM)+1,nrow = 0))
pred_score<-as.data.frame(matrix(ncol=nrow(test_MM)+1,nrow = 0))


svmfit<-readRDS("./set_svmfit.rds")
rffit<-readRDS("./set_rffit.rds")
nn_fit<-readRDS("./set_nn_fit.rds")
lvq_fit<-readRDS("./set_lvq_fit.rds")
gbm_fit<-readRDS("./set_gbm_fit.rds")
glmfit<-readRDS("./set_glmfit.rds")
fit_stacked<-readRDS("./set_fit_stacked.rds")


##SVM
##the random number seed is set before each algorithm is trained to ensure that each algorithm gets the same data partitions and repeats
library(e1071)
tc<-tune.control(cross=10)
svmtrain<-train_MM
svmfit<-svm(x=svmtrain[,1:nfea],y=as.factor(svmtrain$set),gamma=0.001,tune.control=tc,kernel = 'radial',probability=TRUE)
svm_pred_train<-predict(svmfit,svmtrain[,1:nfea])
svm_pred_train<-as.matrix(svm_pred_train)
rownames(svm_pred_train)<-rownames(train_MM)
prob<-predict(svmfit,test_MM[,1:nfea],probability=TRUE)
pred<-predict(svmfit,test_MM[,1:nfea])
test_MM$set<-as.vector(test_MM$set)
pred<-as.matrix(pred)
rownames(pred)<-rownames(test_MM)
SVM_accuracy<-Accuracy(pred,test_MM$set)
SVM_sensitivity<-Sensitivity(test_MM$set,pred)
SVM_classloss<-ZeroOneLoss(test_MM$set,pred)
#add model svm prediction score to class_result and ae_predict
asvm<-as.data.frame(t(c(SVM_accuracy,SVM_classloss)))
asvm<-as.data.frame(round(asvm,4))
asvm.model<-paste("SVM",sep="")
asvm<-cbind(asvm.model,asvm)
colnames(asvm)<-colnames(ae_predict)
ae_predict<-rbind(ae_predict,asvm)
pred.model<-paste("SVM",sep="")
pred_class<-cbind(pred.model,t(pred))
colnames(pred_class)<-colnames(class_result)
class_result<-rbind(class_result,pred_class)


##random forest
set.seed(7)
library(randomForest)
fitcontrol<-trainControl(method="cv",number=10,classProbs = TRUE)
RFtrain<-train_MM
RFtrain$set<-as.factor(RFtrain$set)
levels(RFtrain$set)<-make.names(levels(factor(RFtrain$set)))
rffit<-caret::train(x=RFtrain[,1:nfea],
                    y=RFtrain$set,
                    method="parRF",
                    tuneGrid=expand.grid(
                      .mtry=c(5,2,1)),
                    metric="Accuracy",
                    exclude=NULL,
                    trControl=fitcontrol)
RFtest<-test_MM
RFtest$set<-as.factor(RFtest$set)
levels(RFtest$set)<-make.names(levels(factor(RFtest$set)))
prob_rf<-predict(rffit,RFtest[,1:nfea],type='prob')
pred_rf<-array(1:nrow(prob_rf))
for ( i in 1:nrow(prob_rf)) {
  pred_rf[i]<-as.character(names(which.max(prob_rf[i,]))[1])
} 
pred_rf<-as.matrix(pred_rf)
rownames(pred_rf)<-rownames(test_MM)
rf_prob_train<-predict(rffit,RFtrain[,1:nfea],type="prob")
rf_pred_train<-array(1:nrow(rf_prob_train))
for ( i in 1:nrow(rf_prob_train)) {
  rf_pred_train[i]<-as.character(names(which.max(rf_prob_train[i,]))[1])
}
rf_pred_train<-as.matrix(rf_pred_train)
rownames(rf_pred_train)<-rownames(train_MM)
rf_accuracy<-Accuracy(pred_rf,RFtest$set)
rf_classloss<-ZeroOneLoss(RFtest$set,pred_rf)
#add model rf prediction score to class_result and ae_predict
arf<-as.data.frame(t(c(rf_accuracy,rf_classloss)))
arf<-as.data.frame(round(arf,4))
arf.model<-paste("RF",sep="")
arf<-cbind(arf.model,arf)
colnames(arf)<-colnames(ae_predict)
ae_predict<-rbind(ae_predict,arf)
pred_rf.model<-paste("RF",sep="")
pred_rf_class<-cbind(pred_rf.model,t(pred_rf))
colnames(pred_rf_class)<-colnames(class_result)
class_result<-rbind(class_result,pred_rf_class)

#Neural Network
nncontrol=trainControl(method = "cv", number = 10) 
NNtrain<-train_MM
nn_Grid<-expand.grid(
  .size=c(50,5),
  .decay=0.00147)
nn_fit<-caret::train(set~.,
                     NNtrain,
                     method="nnet",
                     metric="Accuracy",
                     tuneGrid=nn_Grid,
                     MaxNWts=50000,
                     maxit=100,
                     trace=FALSE)
pred_nn<-predict(nn_fit,test_MM[,1:nfea])
pred_nn<-as.matrix(pred_nn)
rownames(pred_nn)<-rownames(test_MM)
nn_pred_train<-predict(nn_fit,RFtrain[,1:nfea])
nn_pred_train<-as.matrix(nn_pred_train)
rownames(nn_pred_train)<-rownames(train_MM)
prob_nn<-predict(nn_fit,test_MM,type='prob')
nn_prob_train<-predict(nn_fit,RFtrain[,1:nfea],type="prob")
#prob_nn<-predict(nn_fit,test_MM,type='prob')
nn_accuracy<-Accuracy(pred_nn,test_MM$set)
nn_classloss<-ZeroOneLoss(test_MM$set,pred_nn)
#add model rf prediction score to class_result and ae_predict
a_nn<-as.data.frame(t(c(nn_accuracy,nn_classloss)))
a_nn<-as.data.frame(round(a_nn,4))
a_nn.model<-paste("ANN",sep="")
a_nn<-cbind(a_nn.model,a_nn)
colnames(a_nn)<-colnames(ae_predict)
ae_predict<-rbind(ae_predict,a_nn)
pred_nn.model<-paste("ANN",sep="")
pred_nn_class<-cbind(pred_nn.model,t(pred_nn))
colnames(pred_nn_class)<-colnames(class_result)
class_result<-rbind(class_result,pred_nn_class)

#Learning Vector Quantization
set.seed(7)
lvqcontrol<-trainControl(method = "cv", number = 10)
lvq_grid <- expand.grid(size=c(5,10,20,50), k=c(2,4,5,7,10))
LVQtrain<-train_MM
lvq_fit<-caret::train(set~.,
                      LVQtrain,
                      method="lvq",
                      trControl=lvqcontrol,
                      tuneGrid=lvq_grid)
pred_lvq<-predict(lvq_fit,test_MM[,1:nfea])
pred_lvq<-as.matrix(pred_lvq)
rownames(pred_lvq)<-rownames(test_MM)
lvq_pred_train<-predict(lvq_fit,RFtrain[,1:nfea])
lvq_pred_train<-as.matrix(lvq_pred_train)
rownames(lvq_pred_train)<-rownames(train_MM)
lvq_accuracy<-Accuracy(pred_lvq,test_MM$set)
lvq_classloss<-ZeroOneLoss(test_MM$set,pred_lvq)
#add model lvq prediction score to class_result and ae_predict
a_lvq<-as.data.frame(t(c(lvq_accuracy,lvq_classloss)))
a_lvq<-as.data.frame(round(a_lvq,4))
a_lvq.model<-paste("lvq",sep="")
a_lvq<-cbind(a_lvq.model,a_lvq)
colnames(a_lvq)<-colnames(ae_predict)
ae_predict<-rbind(ae_predict,a_lvq)
pred_lvq.model<-paste("lvq",sep="")
pred_lvq_class<-cbind(pred_lvq.model,t(pred_lvq))
colnames(pred_lvq_class)<-colnames(class_result)
class_result<-rbind(class_result,pred_lvq_class)


##GBM
set.seed(7)
gbm_fit<-caret::train(set~.,
                      RFtrain,
                      method="gbm",
                      verbose=FALSE,
                      trControl=fitcontrol)
prob_gbm<-predict(gbm_fit,RFtest[,1:nfea],type='prob')
pred_gbm<-array(1:nrow(prob_gbm))
for ( i in 1:nrow(prob_gbm)) {
  pred_gbm[i]<-as.character(names(which.max(prob_gbm[i,]))[1])
} 
pred_gbm<-as.matrix(pred_gbm)
rownames(pred_gbm)<-rownames(RFtest)

gbm_prob_train<-predict(gbm_fit,RFtrain[,1:nfea],type="prob")
gbm_pred_train<-array(1:nrow(gbm_prob_train))
for ( i in 1:nrow(gbm_prob_train)) {
  gbm_pred_train[i]<-as.character(names(which.max(gbm_prob_train[i,]))[1])
}
gbm_pred_train<-as.matrix(gbm_pred_train)
rownames(gbm_pred_train)<-rownames(train_MM)
gbm_accuracy<-Accuracy(pred_gbm,RFtest$set)
gbm_classloss<-ZeroOneLoss(RFtest$set,pred_gbm)
#add model gbm prediction score to class_result and ae_predict
a_gbm<-as.data.frame(t(c(gbm_accuracy,gbm_classloss)))
a_gbm<-as.data.frame(round(a_gbm,4))
a_gbm.model<-paste("GBM",sep="")
a_gbm<-cbind(a_gbm.model,a_gbm)
colnames(a_gbm)<-colnames(ae_predict)
ae_predict<-rbind(ae_predict,a_gbm)
pred_gbm.model<-paste("GBM",sep="")
pred_gbm_class<-cbind(pred_gbm.model,t(pred_gbm))
colnames(pred_gbm_class)<-colnames(class_result)
class_result<-rbind(class_result,pred_gbm_class)


##glmnet Generalized Linear Model
set.seed(7)
gridsearch_for_lambda =  data.frame (alpha = 0,
                                     lambda = c (2^c(-15:15), 3^c(-15:15)))
train_control = trainControl (method="cv", number =10,
                              savePredictions =TRUE , allowParallel = FALSE,classProbs = TRUE )
glmfit<-caret::train(set~.,
                     RFtrain,
                     method="glmnet",
                     tuneGrid=gridsearch_for_lambda,
                     trControl=train_control,
                     preProcess=NULL)
prob_glm<-predict(glmfit,RFtest[,1:nfea],type='prob')
pred_glm<-array(1:nrow(prob_glm))
for ( i in 1:nrow(prob_gbm)) {
  pred_glm[i]<-as.character(names(which.max(prob_glm[i,]))[1])
} 
pred_glm<-as.matrix(pred_glm)
rownames(pred_glm)<-rownames(test_MM)
glm_prob_train<-predict(glmfit,RFtrain[,1:nfea],type="prob")
glm_pred_train<-array(1:nrow(glm_prob_train))
for ( i in 1:nrow(glm_prob_train)) {
  glm_pred_train[i]<-as.character(names(which.max(glm_prob_train[i,]))[1])
} 
glm_pred_train<-as.matrix(glm_pred_train)
rownames(glm_pred_train)<-rownames(train_MM)
glm_accuracy<-Accuracy(pred_glm,RFtest$set)
glm_classloss<-ZeroOneLoss(RFtest$set,pred_glm)
#add model glm prediction score to class_result and ae_predict
a_glm<-as.data.frame(t(c(glm_accuracy,glm_classloss)))
a_glm<-as.data.frame(round(a_glm,4))
a_glm.model<-paste("GLM",sep="")
a_glm<-cbind(a_glm.model,a_glm)
colnames(a_glm)<-colnames(ae_predict)
ae_predict<-rbind(ae_predict,a_glm)
pred_glm.model<-paste("GLM",sep="")
pred_glm_class<-cbind(pred_glm.model,t(pred_glm))
colnames(pred_glm_class)<-colnames(class_result)
class_result<-rbind(class_result,pred_glm_class)

##Stacked model
#combine all the predictions of train data from above classifiers to train stacked model
stacked_tc<-trainControl(method = 'cv',number=10,classProbs = TRUE)
train_set<-RFtrain$set
combo_train<-data.frame(svm_pred_train,rf_pred_train,nn_pred_train,lvq_pred_train,gbm_pred_train,glm_pred_train,train_set)
levels(combo_train$train_set)<-make.names(levels(factor(combo_train$train_set)))
set<-as.factor(test_MM$set)
combo_test<-data.frame(pred,pred_rf,pred_nn,pred_lvq,pred_gbm,pred_glm,set)
colnames(combo_test)<-colnames(combo_train)
levels(combo_test$train_set)<-make.names(levels(factor(combo_test$train_set)))
fit_stacked<-caret::train(as.factor(train_set)~.,
                          combo_train,
                          method="parRF",
                          tuneGrid=expand.grid(
                            .mtry=1),
                          metric="Accuracy",
                          trControl=stacked_tc)
prob_stacked<-predict(fit_stacked,combo_test[,1:6],type='prob')
pred_stacked<-array(1:nrow(prob_stacked))
for ( i in 1:nrow(prob_stacked)) {
  pred_stacked[i]<-as.character(names(which.max(prob_stacked[i,]))[1])
} 
pred_stacked<-as.matrix(pred_stacked)
rownames(pred_stacked)<-rownames(test_MM)
test_MM$set<-as.vector(test_MM$set)
stacked_accuracy<-Accuracy(pred_stacked,combo_test$train_set)
stacked_classloss<-ZeroOneLoss(combo_test$train_set,pred_stacked)
#add model rf prediction score to class_result and ae_predict
a_stacked<-as.data.frame(t(c(stacked_accuracy,stacked_classloss)))
a_stacked<-as.data.frame(round(a_stacked,4))
a_stacked.model<-paste("stacked_model",sep="")
a_stacked<-cbind(a_stacked.model,a_stacked)
colnames(a_stacked)<-colnames(ae_predict)
ae_predict<-rbind(ae_predict,a_stacked)
pred_stacked.model<-paste("stacked_model",sep="")
pred_stacked_class<-cbind(pred_stacked.model,t(pred_stacked))
colnames(pred_stacked_class)<-colnames(class_result)
class_result<-rbind(class_result,pred_stacked_class)

#save scale parameters of set predictive mdoel
write.csv(set_mean,"./set_mean.csv")
write.csv(set_sd,"./set_sd.csv")


#save model seperately
saveRDS(svmfit,"./set_svmfit.rds")
saveRDS(rffit,"./set_rffit.rds")
saveRDS(nn_fit,"./set_nn_fit.rds")
saveRDS(lvq_fit,"./set_lvq_fit.rds")
saveRDS(gbm_fit,"./set_gbm_fit.rds")
saveRDS(glmfit,"./set_glmfit.rds")
saveRDS(fit_stacked,"./set_fit_stacked.rds")


set_svmfit<-readRDS("./set_svmfit.rds")
set_rffit<-readRDS("./set_rffit.rds")
set_nn_fit<-readRDS("./set_nn_fit.rds")
set_lvq_fit<-readRDS("./set_lvq_fit.rds")
set_gbm_fit<-readRDS("./set_gbm_fit.rds")
set_glmfit<-readRDS("./set_glmfit.rds")
set_fit_stacked<-readRDS("./set_fit_stacked.rds")


train_MM<-train_data
test_MM<-test_data
train_MM_set<-train_MM$set
train_MM$set<-NULL
test_MM_set<-test_MM$set
test_MM$set<-NULL



#H_mean<-readRDS("./H_mean.rds")
#H_sd<-readRDS("./H_sd.rds")
M_mean<-readRDS("./M_mean.rds")
M_sd<-readRDS("./M_sd.rds")
E_mean<-readRDS("./E_mean.rds")
E_sd<-readRDS("./E_sd.rds")
U_mean<-readRDS("./U_mean.rds")
U_sd<-readRDS("./U_sd.rds")


## Read in the trained model state.
U_svmfit<-readRDS("./U_svmfit.rds")
U_rffit<-readRDS("./U_rffit.rds")
U_nn_fit<-readRDS("./U_nn_fit.rds")
U_lvq_fit<-readRDS("./U_lvq_fit.rds")
U_gbm_fit<-readRDS("./U_gbm_fit.rds")
U_glmfit<-readRDS("./U_glmfit.rds")
U_fit_stacked<-readRDS("./U_fit_stacked.rds")

E_svmfit<-readRDS("./E_svmfit.rds")
E_rffit<-readRDS("./E_rffit.rds")
E_nn_fit<-readRDS("./E_nn_fit.rds")
E_lvq_fit<-readRDS("./E_lvq_fit.rds")
E_gbm_fit<-readRDS("./E_gbm_fit.rds")
E_glmfit<-readRDS("./E_glmfit.rds")
E_fit_stacked<-readRDS("./E_fit_stacked.rds")

#H_svmfit<-readRDS("./H_svmfit.rds")
#H_rffit<-readRDS("./H_rffit.rds")
#H_nn_fit<-readRDS("./H_nn_fit.rds")
#H_lvq_fit<-readRDS("./H_lvq_fit.rds")
#H_gbm_fit<-readRDS("./H_gbm_fit.rds")
#H_glmfit<-readRDS("./H_glmfit.rds")
#H_fit_stacked<-readRDS("./H_fit_stacked.rds")

M_svmfit<-readRDS("./M_svmfit.rds")
M_rffit<-readRDS("./M_rffit.rds")
M_nn_fit<-readRDS("./M_nn_fit.rds")
M_lvq_fit<-readRDS("./M_lvq_fit.rds")
M_gbm_fit<-readRDS("./M_gbm_fit.rds")
M_glmfit<-readRDS("./M_glmfit.rds")
M_fit_stacked<-readRDS("./M_fit_stacked.rds")

#eu_svmfit<-readRDS("./eu_svmfit.rds")
#eu_rffit<-readRDS("./eu_rffit.rds")
#eu_nn_fit<-readRDS("./eu_nn_fit.rds")
#eu_lvq_fit<-readRDS("./eu_lvq_fit.rds")
#eu_gbm_fit<-readRDS("./eu_gbm_fit.rds")
#eu_glmfit<-readRDS("./eu_glmfit.rds")
#eu_fit_stacked<-readRDS("./eu_fit_stacked.rds")

risk_flag<-setNames(as.data.frame(matrix(0,ncol=1,nrow = nrow(test_MM))),c("risk_flag"))
rownames(risk_flag)<-rownames(test_MM)
th<-0.5

##according to the predict set each sample belong to, feed into sepecifc dataset model
for (i in  1: nrow(pred_stacked)) {
  if (pred_stacked[i]=="X1") {
      input<-test_MM[i,]
      test_U<-setNames(as.data.frame(matrix(0,nrow=1,ncol=length(U_genes_eu_list))),U_genes_eu_list)
      idx<-colnames(input)[colnames(input)%in%colnames(test_U)]
      test_U[,idx]<-input[,idx]
      test=sweep(test_U,MARGIN=2,U_mean,"-")
      test_U<-sweep(test,MARGIN=2,U_sd,"/")
      test_U$label<-input$label
      ##feed into U_models
      test_U_svm_prob<-predict(U_svmfit,test_U[,1:(ncol(test_U)-1)],probability = TRUE)
      test_U_svm<-predict(U_svmfit,test_U[,1:(ncol(test_U)-1)])
      test_U_rf_prob<-predict(U_rffit,test_U[,1:(ncol(test_U)-1)],type="prob")
      test_U_rf<-factor(ifelse(test_U_rf_prob$X1>th,"1","-1"))
      test_U_nn_prob<-predict(U_nn_fit,test_U[,1:(ncol(test_U)-1)],type="prob")
      test_U_nn<-predict(U_nn_fit,test_U[,1:(ncol(test_U)-1)])
      test_U_lvq<-predict(U_lvq_fit,test_U[,1:(ncol(test_U)-1)])
      test_U_gbm_prob<-predict(U_gbm_fit,test_U[,1:(ncol(test_U)-1)],type="prob")
      test_U_gbm<-factor(ifelse(test_U_gbm_prob$X1>th,"1","-1"))
      test_U_glm_prob<-predict(U_glmfit,test_U[,1:(ncol(test_U)-1)],type="prob")
      test_U_glm<-factor(ifelse(test_U_glm_prob$X1>th,"1","-1"))
      combo_U_test<-data.frame(test_U_svm,test_U_rf,test_U_nn,test_U_lvq,test_U_gbm,test_U_glm,label<-as.factor(test_U$label))
      colnames(combo_U_test)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train","train_label" )
      test_U_prob<-predict(U_fit_stacked,combo_U_test,type="prob")
      test_U_stacked<-factor(ifelse(test_U_prob$X1>th,"1","-1"))
      if (test_U_stacked=="1"){test_U_stacked=1} else {test_U_stacked=-1}
      risk_flag[i,1]<-test_U_stacked
  }
  else if (pred_stacked[i]=="X2") {
    input<-test_MM[i,]
    test_E<-setNames(as.data.frame(matrix(0,nrow=1,ncol=length(E_genes_eu_list))),E_genes_eu_list)
    idx<-colnames(input)[colnames(input)%in%colnames(test_E)]
    test_E[,idx]<-input[,idx]
    test=sweep(test_E,MARGIN=2,E_mean,"-")
    test_E<-sweep(test,MARGIN=2,E_sd,"/")
    test_E$label<-input$label
    ##load the E_mdoel
    test_E_svm_prob<-predict(E_svmfit,test_E[,1:(ncol(test_E)-1)],probability = TRUE)
    test_E_svm<-predict(E_svmfit,test_E[,1:(ncol(test_E)-1)])
    test_E_rf_prob<-predict(E_rffit,test_E[,1:(ncol(test_E)-1)],type="prob")
    test_E_rf<-factor(ifelse(test_E_rf_prob$X1>th,"1","-1"))
    test_E_nn_prob<-predict(E_nn_fit,test_E[,1:(ncol(test_E)-1)],type="prob")
    test_E_nn<-predict(E_nn_fit,test_E[,1:(ncol(test_E)-1)])
    test_E_lvq<-predict(E_lvq_fit,test_E[,1:(ncol(test_E)-1)])
    test_E_gbm_prob<-predict(E_gbm_fit,test_E[,1:(ncol(test_E)-1)],type="prob")
    test_E_gbm<-factor(ifelse(test_E_gbm_prob$X1>th,"1","-1"))
    test_E_glm_prob<-predict(E_glmfit,test_E[,1:(ncol(test_E)-1)],type="prob")
    test_E_glm<-factor(ifelse(test_E_glm_prob$X1>th,"1","-1"))
    combo_E_test<-data.frame(test_E_svm,test_E_rf,test_E_nn,test_E_lvq,test_E_gbm,test_E_glm,label<-as.factor(test_E$label))
    colnames(combo_E_test)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train" ,"train_label" )
    test_E_prob<-predict(E_fit_stacked,combo_E_test,type="prob")
    th<-0.5
    test_E_stacked<-factor(ifelse(test_E_prob$X1>th,"1","-1"))
    if (test_E_stacked=="1"){test_E_stacked=1} else {test_E_stacked=-1}
    risk_flag[i,1]<-test_E_stacked
  }
  else if (pred_stacked[i]=="X3") {
    input<-test_MM[i,]
    test_M<-setNames(as.data.frame(matrix(0,nrow=1,ncol=length(M_genes_eu_list))),M_genes_eu_list)
    idx<-colnames(input)[colnames(input)%in%colnames(test_M)]
    test_M[,idx]<-input[,idx]
    test=sweep(test_M,MARGIN=2,M_mean,"-")
    test_M<-sweep(test,MARGIN=2,M_sd,"/")
    test_M$label<-input$label
    ##load the M_mdoel
    test_M_svm<-predict(M_svmfit,test_M[,1:(ncol(test_M)-1)])
    test_M_rf_prob<-predict(M_rffit,test_M[,1:(ncol(test_M)-1)],type="prob")
    test_M_rf<-factor(ifelse(test_M_rf_prob$X1>th,"1","-1"))
    test_M_nn_prob<-predict(M_nn_fit,test_M[,1:(ncol(test_M)-1)],type="prob")
    test_M_nn<-predict(M_nn_fit,test_M[,1:(ncol(test_M)-1)])
    test_M_lvq<-predict(M_lvq_fit,test_M[,1:(ncol(test_M)-1)])
    test_M_gbm_prob<-predict(M_gbm_fit,test_M[,1:(ncol(test_M)-1)],type="prob")
    test_M_gbm<-factor(ifelse(test_M_gbm_prob$X1>th,"1","-1"))
    test_M_glm_prob<-predict(M_glmfit,test_M[,1:(ncol(test_M)-1)],type="prob")
    test_M_glm<-factor(ifelse(test_M_glm_prob$X1>th,"1","-1"))
    combo_M_test<-data.frame(test_M_svm,test_M_rf,test_M_nn,test_M_lvq,test_M_gbm,test_M_glm,label<-as.factor(test_M$label))
    colnames(combo_M_test)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train","train_label" )
    test_M_prob<-predict(M_fit_stacked,combo_M_test,type="prob")
    test_M_stacked<-factor(ifelse(test_M_prob$X1>th,"1","-1"))
    if (test_M_stacked=="1"){test_M_stacked=1} else {test_M_stacked=-1}
    risk_flag[i,1]<-test_M_stacked
  }
}


# Area under the curve: 0.935


save(U_svmfit,U_rffit,U_nn_fit,U_lvq_fit,U_gbm_fit,U_glmfit,U_fit_stacked,
     E_svmfit,E_rffit,E_nn_fit,E_lvq_fit,E_gbm_fit,E_glmfit,E_fit_stacked,
     M_svmfit,M_rffit,M_nn_fit,M_lvq_fit,M_gbm_fit,M_glmfit,M_fit_stacked,
     set_svmfit,set_rffit,set_nn_fit,set_lvq_fit,set_gbm_fit,set_glmfit,set_fit_stacked,
     U_genes_eu_list,E_genes_eu_list,M_genes_eu_list,
     U_genes_list,E_genes_list,H_genes_list,M_genes_list,eu_genes_list,
     set_mean,set_sd,M_mean,M_sd,U_mean,U_sd,E_mean,E_sd,
     file="model-state-metadata.Rd")

