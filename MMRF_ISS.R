library(gbm)
library(MLmetrics)
library(ModelMetrics)
library(caret)
library(plyr)
library(FCNN4R)
library(dplyr)



setwd("/Users/guanhaibin/MMRF_ISS")
M_genes_list<-read.csv("M_genes_eu.csv",row.names = 1)
M_genes_list<-as.character(M_genes_list[,1])
model.entrezs<-M_genes_list

#read MMRF datasets with selected genes
Gene_expression=read.csv("MMRFentrezIDlevel.csv",header=T,sep=",",row.names=1)
Gene_expression<-as.data.frame(Gene_expression)
colnames(Gene_expression)<-gsub("X","",colnames(Gene_expression))

clinical<-read.csv("globalClinTraining.csv",header=T,sep=",",row.names = 2,stringsAsFactors = FALSE)
a<-c("MMRF")
clinical<-subset(clinical,Study%in%a)
clinical<-subset(clinical[,c("D_Age","D_ISS","HR_FLAG")])
clinical[clinical[,"HR_FLAG"]==TRUE,3]=1
clinical[clinical[,"HR_FLAG"]==FALSE,3]=-1
gene_clinical<-subset(clinical,rownames(clinical)%in%rownames(Gene_expression))

select.data<-setNames(as.data.frame(matrix(0,nrow=nrow(Gene_expression),ncol=length(model.entrezs))),model.entrezs)
idx<-colnames(Gene_expression)[colnames(Gene_expression)%in%colnames(select.data)]
select.data[,idx]<-Gene_expression[,idx]
rownames(select.data)<-rownames(Gene_expression)
expression <- select.data

train_data<-merge(expression,gene_clinical,by="row.names",all.x=TRUE)
row.names(train_data)<-train_data[,1]
train_data<-train_data[-c(1)]
colnames(train_data)[ncol(train_data)]<-c("label")
train_data$label<-as.factor(train_data$label)
train_data<-subset(train_data,D_ISS!="NA")
train_data$D_Age<-train_data$D_Age/10


####ISS=I
train_ISS_I<-subset(train_data,D_ISS=="1")
b<-c("D_ISS")
train_ISS_I<-train_ISS_I[,!(names(train_ISS_I)%in%b)]
data_I<-ubBalance(X=train_ISS_I[,1:(ncol(train_ISS_I)-1)], Y=train_ISS_I$label, type="ubSMOTE", positive=1, percOver=100, percUnder=200,
                k=5, perc=50, method="percPos", w=NULL, verbose=FALSE)
data_I<-cbind(data_I$X,data_I$Y)
train_I<-data_I
colnames(train_I)[ncol(train_I)]<-c("label")
set.seed(123)
train.rows<-createDataPartition(y=train_I$label,p=0.85,list=FALSE)
train_MM<-train_I[train.rows,]
test_MM<-train_I[-train.rows,]

nfea<-ncol(train_MM)-1




##create a dataframe to store prediction score for different model
ae_predict<-setNames(as.data.frame(matrix(ncol=9,nrow = 0)),c("model","Accuracy","AUC","F1_score","PRAUC","MCC","sensitivity","specific","classloss"))
class_result<-as.data.frame(matrix(ncol=nrow(test_MM)+1,nrow = 0))
pred_score<-as.data.frame(matrix(ncol=nrow(test_MM)+1,nrow = 0))
#function change label to 0 and 1 for calculate AUC and MMC
adjustAUC<-function(predict,truth){
  {for (i in 1:length(predict)) {if (predict[i]==-1 ){{predict[i]=0}}}}
  {for (i in 1:length(truth)) {if (truth[i]==-1 ){{truth[i]=0}}}}
  return(AUC(predict,truth))
}

adjustmcc<-function(truth,predict,cutoff=1){
  {for (i in 1:length(predict)) {if (predict[i]==-1 ){{predict[i]=0}}}}
  {for (i in 1:length(truth)) {if (truth[i]==-1 ){{truth[i]=0}}}}
  return(mcc(truth,predict,cutoff=1))
}


##SVM
##the random number seed is set before each algorithm is trained to ensure that each algorithm gets the same data partitions and repeats
library(e1071)
tc<-tune.control(cross=10)
svmtrain<-train_MM
svmfit<-svm(x=svmtrain[,1:nfea],y=as.factor(svmtrain$label),gamma=0.01,tune.control=tc,scale=TRUE,kernel = 'radial',probability=TRUE)
svm_pred_train<-predict(svmfit,svmtrain[,1:nfea])
prob<-predict(svmfit,test_MM[,1:nfea],probability=TRUE)
pred<-predict(svmfit,test_MM[,1:nfea])
pred<-as.vector(pred)
test_MM$label<-as.vector(test_MM$label)
SVM_accuracy<-Accuracy(pred,test_MM$label)
SVM_AUC<-adjustAUC(as.numeric(pred),as.numeric(test_MM$label))
SVM_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred))
SVM_PRAUC<-PRAUC(as.numeric(pred),as.numeric(test_MM$label))
SVM_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred),cutoff=1)
SVM_sensitivity<-Sensitivity(test_MM$label,pred)
SVM_specific<-Specificity(test_MM$label,pred)
SVM_classloss<-ZeroOneLoss(test_MM$label,pred)
#add model svm prediction score to class_result and ae_predict
asvm<-as.data.frame(t(c(SVM_accuracy,SVM_AUC,SVM_F1,SVM_PRAUC,SVM_MCC,SVM_sensitivity,SVM_specific,SVM_classloss)))
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
levels(RFtrain$label)<-make.names(levels(factor(RFtrain$label)))
rffit<-caret::train(x=RFtrain[,1:nfea],
                    y=RFtrain$label,
                    method="parRF",
                    tuneGrid=expand.grid(
                      .mtry=c(5,3,2,1)),
                    metric="Accuracy",
                    preProcess = c('center', 'scale'),
                    exclude=NULL,
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
rf_AUC<-adjustAUC(as.numeric(pred_rf),as.numeric(test_MM$label))
rf_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_rf))
rf_PRAUC<-PRAUC(as.numeric(pred_rf),as.numeric(test_MM$label))
rf_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_rf),cutoff=1)
rf_sensitivity<-Sensitivity(test_MM$label,pred_rf)
rf_specific<-Specificity(test_MM$label,pred_rf)
rf_classloss<-ZeroOneLoss(test_MM$label,pred_rf)
#add model rf prediction score to class_result and ae_predict
arf<-as.data.frame(t(c(rf_accuracy,rf_AUC,rf_F1,rf_PRAUC,rf_MCC,rf_sensitivity,rf_specific,rf_classloss)))
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
nn_fit<-caret::train(label~.,
                     NNtrain,
                     method="nnet",
                     metric="Accuracy",
                     preProcess = c('center', 'scale'),
                     tuneGrid=nn_Grid,
                     MaxNWts=10000,
                     maxit=100,
                     trace=FALSE)
pred_nn<-predict(nn_fit,test_MM[,1:nfea])
nn_pred_train<-predict(nn_fit,RFtrain[,1:nfea])
prob_nn<-predict(nn_fit,test_MM,type='prob')
nn_prob_train<-predict(nn_fit,RFtrain[,1:nfea],type="prob")
nn_pred_train<-predict(nn_fit,RFtrain[,1:nfea])
#prob_nn<-predict(nn_fit,test_MM,type='prob')
pred_nn<-as.vector(pred_nn)
test_MM$label<-as.vector(test_MM$label)
nn_accuracy<-Accuracy(pred_nn,test_MM$label)
nn_AUC<-adjustAUC(as.numeric(pred_nn),as.numeric(test_MM$label))
nn_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_nn))
nn_PRAUC<-PRAUC(as.numeric(pred_nn),as.numeric(test_MM$label))
nn_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_nn),cutoff=1)
nn_sensitivity<-Sensitivity(test_MM$label,pred_nn)
nn_specific<-Specificity(test_MM$label,pred_nn)
nn_classloss<-ZeroOneLoss(test_MM$label,pred_nn)
#add model rf prediction score to class_result and ae_predict
a_nn<-as.data.frame(t(c(nn_accuracy,nn_AUC,nn_F1,nn_PRAUC,nn_MCC,nn_sensitivity,nn_specific,nn_classloss)))
a_nn<-as.data.frame(round(a_nn,4))
a_nn.model<-paste("ANN",sep="")
a_nn<-cbind(a_nn.model,a_nn)
colnames(a_nn)<-colnames(ae_predict)
ae_predict<-rbind(ae_predict,a_nn)
pred_nn.model<-paste("ANN",sep="")
pred_nn_class<-cbind(pred_nn.model,t(pred_nn))
colnames(pred_nn_class)<-colnames(class_result)
class_result<-rbind(class_result,pred_nn_class)
#Learning Vector Quantization -> xgboost
set.seed(7)
lvqcontrol<-trainControl(method = "cv", number = 3, returnResamp = "all", search = "random",classProbs = TRUE, summaryFunction = twoClassSummary)
LVQtrain<-train_MM
lvq_grid <-expand.grid(nrounds = c(1, 10),
                       max_depth = c(1, 4),
                       eta = c(.1, .4),
                       gamma = 0,
                       colsample_bytree = .7,
                       min_child_weight = 1,
                       subsample = c(.8, 1))
lvq_fit<-caret::train(label~.,
                      RFtrain,
                      method="xgbTree",
                      trControl=lvqcontrol,
                      preProcess = c('center', 'scale'),
                      turnGrid=lvq_grid,
                      metric="ROC",
                      verbose=FALSE)
prob_lvq<-predict(lvq_fit,test_MM[,1:nfea],type='prob')
th<-0.5
pred_lvq<-factor(ifelse(prob_lvq$X1>th,"1","-1"))
lvq_prob_train<-predict(lvq_fit,RFtrain[,1:nfea],type="prob")
lvq_pred_train<-factor(ifelse(lvq_prob_train$X1>th,"1","-1"))
pred_lvq<-as.vector(pred_lvq)
test_MM$label<-as.vector(test_MM$label)
lvq_accuracy<-Accuracy(pred_lvq,test_MM$label)
lvq_AUC<-adjustAUC(pred_lvq,as.numeric(test_MM$label))
lvq_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_lvq))
lvq_PRAUC<-PRAUC(as.numeric(pred_lvq),as.numeric(test_MM$label))
lvq_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_lvq),cutoff=1)
lvq_sensitivity<-Sensitivity(test_MM$label,pred_lvq)
lvq_specific<-Specificity(test_MM$label,pred_lvq)
lvq_classloss<-ZeroOneLoss(test_MM$label,pred_lvq)
#add model lvq prediction score to class_result and ae_predict
a_lvq<-as.data.frame(t(c(lvq_accuracy,lvq_AUC,lvq_F1,lvq_PRAUC,lvq_MCC,lvq_sensitivity,lvq_specific,lvq_classloss)))
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
gbm_fit<-caret::train(label~.,
                      RFtrain,
                      method="gbm",
                      preProcess = c('center', 'scale'),
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
gbm_AUC<-adjustAUC(pred_gbm,as.numeric(test_MM$label))
gbm_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_gbm))
gbm_PRAUC<-PRAUC(as.numeric(pred_gbm),as.numeric(test_MM$label))
gbm_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_gbm),cutoff=1)
gbm_sensitivity<-Sensitivity(test_MM$label,pred_gbm)
gbm_specific<-Specificity(test_MM$label,pred_gbm)
gbm_classloss<-ZeroOneLoss(test_MM$label,pred_gbm)
#add model gbm prediction score to class_result and ae_predict
a_gbm<-as.data.frame(t(c(gbm_accuracy,gbm_AUC,gbm_F1,gbm_PRAUC,gbm_MCC,gbm_sensitivity,gbm_specific,gbm_classloss)))
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
gridsearch_for_lambda =  data.frame (alpha =  seq(.05, 1, length = 15),
                                     lambda = c ((1:5)/10))
train_control = trainControl (method = "LOOCV", returnResamp = "all", search = "random",classProbs = TRUE,summaryFunction = twoClassSummary )
glmfit<-caret::train(label~.,
                     RFtrain,
                     method="glmnet",
                     tuneGrid=gridsearch_for_lambda,
                     preProcess = c('center', 'scale'),
                     metric="ROC",
                     trControl=train_control)
prob_glm<-predict(glmfit,test_MM[,1:nfea],type='prob')
th<-0.5
pred_glm<-factor(ifelse(prob_glm$X1>th,"1","-1"))
glm_prob_train<-predict(glmfit,RFtrain[,1:nfea],type="prob")
glm_pred_train<-factor(ifelse(glm_prob_train$X1>th,"1","-1"))
pred_glm<-as.vector(pred_glm)
test_MM$label<-as.vector(test_MM$label)
glm_accuracy<-Accuracy(pred_glm,test_MM$label)
glm_AUC<-adjustAUC(pred_glm,as.numeric(test_MM$label))
glm_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_glm))
glm_PRAUC<-PRAUC(as.numeric(pred_glm),as.numeric(test_MM$label))
glm_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_glm),cutoff=1)
glm_sensitivity<-Sensitivity(test_MM$label,pred_glm)
glm_specific<-Specificity(test_MM$label,pred_glm)
glm_classloss<-ZeroOneLoss(test_MM$label,pred_glm)
#add model rf prediction score to class_result and ae_predict
a_glm<-as.data.frame(t(c(glm_accuracy,glm_AUC,glm_F1,glm_PRAUC,glm_MCC,glm_sensitivity,glm_specific,glm_classloss)))
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
train_label<-RFtrain$label
combo_train<-data.frame(svm_pred_train,rf_pred_train,nn_pred_train,lvq_pred_train,gbm_pred_train,glm_pred_train,train_label)
levels(combo_train$train_label)<-make.names(levels(factor(combo_train$train_label)))
label<-as.factor(test_MM$label)
combo_test<-data.frame(pred,pred_rf,pred_nn,pred_lvq,pred_gbm,pred_glm,label)
colnames(combo_test)<-colnames(combo_train) 
levels(combo_test$train_label)<-make.names(levels(factor(combo_test$train_label)))
fit_stacked<-caret::train(as.factor(train_label)~.,
                          combo_train,
                          method="parRF",
                          tuneGrid=expand.grid(
                            .mtry=c(3,2,1)),
                          metric="Accuracy",
                          trControl=stacked_tc)
prob_stacked<-predict(fit_stacked,combo_test[,1:7],type='prob')
th<-0.5
pred_stacked<-factor(ifelse(prob_stacked$X1>th,"1","-1"))
pred_stacked<-as.vector(pred_stacked)
test_MM$label<-as.vector(test_MM$label)
stacked_accuracy<-Accuracy(pred_stacked,test_MM$label)
stacked_AUC<-adjustAUC(pred_stacked,as.numeric(test_MM$label))
stacked_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_stacked))
stacked_PRAUC<-PRAUC(as.numeric(pred_stacked),as.numeric(test_MM$label))
stacked_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_stacked),cutoff=1)
stacked_sensitivity<-Sensitivity(test_MM$label,pred_stacked)
stacked_specific<-Specificity(test_MM$label,pred_stacked)
stacked_classloss<-ZeroOneLoss(test_MM$label,pred_stacked)
#add model rf prediction score to class_result and ae_predict
a_stacked<-as.data.frame(t(c(stacked_accuracy,stacked_AUC,stacked_F1,stacked_PRAUC,stacked_MCC,stacked_sensitivity,stacked_specific,stacked_classloss)))
a_stacked<-as.data.frame(round(a_stacked,4))
a_stacked.model<-paste("stacked_model",sep="")
a_stacked<-cbind(a_stacked.model,a_stacked)
colnames(a_stacked)<-colnames(ae_predict)
ae_predict<-rbind(ae_predict,a_stacked)
pred_stacked.model<-paste("stacked_model",sep="")
pred_stacked_class<-cbind(pred_stacked.model,t(pred_stacked))
colnames(pred_stacked_class)<-colnames(class_result)
class_result<-rbind(class_result,pred_stacked_class)

#save model seperately
saveRDS(svmfit,"./I_svmfit.rds")
saveRDS(rffit,"./I_rffit.rds")
saveRDS(nn_fit,"./I_nn_fit.rds")
saveRDS(lvq_fit,"./I_lvq_fit.rds")
saveRDS(gbm_fit,"./I_gbm_fit.rds")
saveRDS(glmfit,"./I_glmfit.rds")
saveRDS(fit_stacked,"./I_fit_stacked.rds")


####ISS=II
train_ISS_II<-subset(train_data,D_ISS=="2")
b<-c("D_ISS")
train_ISS_II<-train_ISS_II[,!(names(train_ISS_II)%in%b)]
data_II<-ubBalance(X=train_ISS_II[,1:(ncol(train_ISS_II)-1)], Y=train_ISS_II$label, type="ubSMOTE", positive=1, percOver=100, percUnder=200,
                  k=5, perc=50, method="percPos", w=NULL, verbose=FALSE)
data_II<-cbind(data_II$X,data_II$Y)
train_II<-data_II
colnames(train_II)[ncol(train_II)]<-c("label")
set.seed(123)
train.rows<-createDataPartition(y=train_II$label,p=0.85,list=FALSE)
train_MM<-train_II[train.rows,]
test_MM<-train_II[-train.rows,]

nfea<-ncol(train_MM)-1





##create a dataframe to store prediction score for different model
ae_predict<-setNames(as.data.frame(matrix(ncol=9,nrow = 0)),c("model","Accuracy","AUC","F1_score","PRAUC","MCC","sensitivity","specific","classloss"))
class_result<-as.data.frame(matrix(ncol=nrow(test_MM)+1,nrow = 0))
pred_score<-as.data.frame(matrix(ncol=nrow(test_MM)+1,nrow = 0))
#function change label to 0 and 1 for calculate AUC and MMC
adjustAUC<-function(predict,truth){
  {for (i in 1:length(predict)) {if (predict[i]==-1 ){{predict[i]=0}}}}
  {for (i in 1:length(truth)) {if (truth[i]==-1 ){{truth[i]=0}}}}
  return(AUC(predict,truth))
}

adjustmcc<-function(truth,predict,cutoff=1){
  {for (i in 1:length(predict)) {if (predict[i]==-1 ){{predict[i]=0}}}}
  {for (i in 1:length(truth)) {if (truth[i]==-1 ){{truth[i]=0}}}}
  return(mcc(truth,predict,cutoff=1))
}


##SVM
##the random number seed is set before each algorithm is trained to ensure that each algorithm gets the same data partitions and repeats
library(e1071)
tc<-tune.control(cross=10)
svmtrain<-train_MM
svmfit<-svm(x=svmtrain[,1:nfea],y=as.factor(svmtrain$label),gamma=0.01,tune.control=tc,scale=TRUE,kernel = 'radial',probability=TRUE)
svm_pred_train<-predict(svmfit,svmtrain[,1:nfea])
prob<-predict(svmfit,test_MM[,1:nfea],probability=TRUE)
pred<-predict(svmfit,test_MM[,1:nfea])
pred<-as.vector(pred)
test_MM$label<-as.vector(test_MM$label)
SVM_accuracy<-Accuracy(pred,test_MM$label)
SVM_AUC<-adjustAUC(as.numeric(pred),as.numeric(test_MM$label))
SVM_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred))
SVM_PRAUC<-PRAUC(as.numeric(pred),as.numeric(test_MM$label))
SVM_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred),cutoff=1)
SVM_sensitivity<-Sensitivity(test_MM$label,pred)
SVM_specific<-Specificity(test_MM$label,pred)
SVM_classloss<-ZeroOneLoss(test_MM$label,pred)
#add model svm prediction score to class_result and ae_predict
asvm<-as.data.frame(t(c(SVM_accuracy,SVM_AUC,SVM_F1,SVM_PRAUC,SVM_MCC,SVM_sensitivity,SVM_specific,SVM_classloss)))
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
levels(RFtrain$label)<-make.names(levels(factor(RFtrain$label)))
rffit<-caret::train(x=RFtrain[,1:nfea],
                    y=RFtrain$label,
                    method="parRF",
                    tuneGrid=expand.grid(
                      .mtry=c(5,3,2,1)),
                    metric="Accuracy",
                    preProcess = c('center', 'scale'),
                    exclude=NULL,
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
rf_AUC<-adjustAUC(as.numeric(pred_rf),as.numeric(test_MM$label))
rf_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_rf))
rf_PRAUC<-PRAUC(as.numeric(pred_rf),as.numeric(test_MM$label))
rf_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_rf),cutoff=1)
rf_sensitivity<-Sensitivity(test_MM$label,pred_rf)
rf_specific<-Specificity(test_MM$label,pred_rf)
rf_classloss<-ZeroOneLoss(test_MM$label,pred_rf)
#add model rf prediction score to class_result and ae_predict
arf<-as.data.frame(t(c(rf_accuracy,rf_AUC,rf_F1,rf_PRAUC,rf_MCC,rf_sensitivity,rf_specific,rf_classloss)))
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
  .size=c(30,5),
  .decay=0.00147)
nn_fit<-caret::train(label~.,
                     NNtrain,
                     method="nnet",
                     metric="Accuracy",
                     preProcess = c('center', 'scale'),
                     tuneGrid=nn_Grid,
                     MaxNWts=10000,
                     maxit=50,
                     trace=FALSE)
pred_nn<-predict(nn_fit,test_MM[,1:nfea])
nn_pred_train<-predict(nn_fit,RFtrain[,1:nfea])
prob_nn<-predict(nn_fit,test_MM,type='prob')
nn_prob_train<-predict(nn_fit,RFtrain[,1:nfea],type="prob")
nn_pred_train<-predict(nn_fit,RFtrain[,1:nfea])
#prob_nn<-predict(nn_fit,test_MM,type='prob')
pred_nn<-as.vector(pred_nn)
test_MM$label<-as.vector(test_MM$label)
nn_accuracy<-Accuracy(pred_nn,test_MM$label)
nn_AUC<-adjustAUC(as.numeric(pred_nn),as.numeric(test_MM$label))
nn_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_nn))
nn_PRAUC<-PRAUC(as.numeric(pred_nn),as.numeric(test_MM$label))
nn_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_nn),cutoff=1)
nn_sensitivity<-Sensitivity(test_MM$label,pred_nn)
nn_specific<-Specificity(test_MM$label,pred_nn)
nn_classloss<-ZeroOneLoss(test_MM$label,pred_nn)
#add model rf prediction score to class_result and ae_predict
a_nn<-as.data.frame(t(c(nn_accuracy,nn_AUC,nn_F1,nn_PRAUC,nn_MCC,nn_sensitivity,nn_specific,nn_classloss)))
a_nn<-as.data.frame(round(a_nn,4))
a_nn.model<-paste("ANN",sep="")
a_nn<-cbind(a_nn.model,a_nn)
colnames(a_nn)<-colnames(ae_predict)
ae_predict<-rbind(ae_predict,a_nn)
pred_nn.model<-paste("ANN",sep="")
pred_nn_class<-cbind(pred_nn.model,t(pred_nn))
colnames(pred_nn_class)<-colnames(class_result)
class_result<-rbind(class_result,pred_nn_class)
#Learning Vector Quantization -> xgboost
set.seed(7)
lvqcontrol<-trainControl(method = "cv", number = 3, returnResamp = "all", search = "random",classProbs = TRUE, summaryFunction = twoClassSummary)
LVQtrain<-train_MM
lvq_grid <-expand.grid(nrounds = c(1, 10),
                       max_depth = c(1, 4),
                       eta = c(.1, .4),
                       gamma = 0,
                       colsample_bytree = .7,
                       min_child_weight = 1,
                       subsample = c(.8, 1))
lvq_fit<-caret::train(label~.,
                      RFtrain,
                      method="xgbTree",
                      trControl=lvqcontrol,
                      preProcess = c('center', 'scale'),
                      turnGrid=lvq_grid,
                      metric="ROC",
                      verbose=FALSE)
prob_lvq<-predict(lvq_fit,test_MM[,1:nfea],type='prob')
th<-0.5
pred_lvq<-factor(ifelse(prob_lvq$X1>th,"1","-1"))
lvq_prob_train<-predict(lvq_fit,RFtrain[,1:nfea],type="prob")
lvq_pred_train<-factor(ifelse(lvq_prob_train$X1>th,"1","-1"))
pred_lvq<-as.vector(pred_lvq)
test_MM$label<-as.vector(test_MM$label)
lvq_accuracy<-Accuracy(pred_lvq,test_MM$label)
lvq_AUC<-adjustAUC(pred_lvq,as.numeric(test_MM$label))
lvq_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_lvq))
lvq_PRAUC<-PRAUC(as.numeric(pred_lvq),as.numeric(test_MM$label))
lvq_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_lvq),cutoff=1)
lvq_sensitivity<-Sensitivity(test_MM$label,pred_lvq)
lvq_specific<-Specificity(test_MM$label,pred_lvq)
lvq_classloss<-ZeroOneLoss(test_MM$label,pred_lvq)
#add model lvq prediction score to class_result and ae_predict
a_lvq<-as.data.frame(t(c(lvq_accuracy,lvq_AUC,lvq_F1,lvq_PRAUC,lvq_MCC,lvq_sensitivity,lvq_specific,lvq_classloss)))
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
gbm_fit<-caret::train(label~.,
                      RFtrain,
                      method="gbm",
                      preProcess = c('center', 'scale'),
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
gbm_AUC<-adjustAUC(pred_gbm,as.numeric(test_MM$label))
gbm_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_gbm))
gbm_PRAUC<-PRAUC(as.numeric(pred_gbm),as.numeric(test_MM$label))
gbm_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_gbm),cutoff=1)
gbm_sensitivity<-Sensitivity(test_MM$label,pred_gbm)
gbm_specific<-Specificity(test_MM$label,pred_gbm)
gbm_classloss<-ZeroOneLoss(test_MM$label,pred_gbm)
#add model gbm prediction score to class_result and ae_predict
a_gbm<-as.data.frame(t(c(gbm_accuracy,gbm_AUC,gbm_F1,gbm_PRAUC,gbm_MCC,gbm_sensitivity,gbm_specific,gbm_classloss)))
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
gridsearch_for_lambda =  data.frame (alpha =  seq(.05, 1, length = 15),
                                     lambda = c ((1:5)/10))
train_control = trainControl (method = "LOOCV", returnResamp = "all", search = "random",classProbs = TRUE,summaryFunction = twoClassSummary )
glmfit<-caret::train(label~.,
                     RFtrain,
                     method="glmnet",
                     tuneGrid=gridsearch_for_lambda,
                     preProcess = c('center', 'scale'),
                     metric="ROC",
                     trControl=train_control)
prob_glm<-predict(glmfit,test_MM[,1:nfea],type='prob')
th<-0.5
pred_glm<-factor(ifelse(prob_glm$X1>th,"1","-1"))
glm_prob_train<-predict(glmfit,RFtrain[,1:nfea],type="prob")
glm_pred_train<-factor(ifelse(glm_prob_train$X1>th,"1","-1"))
pred_glm<-as.vector(pred_glm)
test_MM$label<-as.vector(test_MM$label)
glm_accuracy<-Accuracy(pred_glm,test_MM$label)
glm_AUC<-adjustAUC(pred_glm,as.numeric(test_MM$label))
glm_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_glm))
glm_PRAUC<-PRAUC(as.numeric(pred_glm),as.numeric(test_MM$label))
glm_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_glm),cutoff=1)
glm_sensitivity<-Sensitivity(test_MM$label,pred_glm)
glm_specific<-Specificity(test_MM$label,pred_glm)
glm_classloss<-ZeroOneLoss(test_MM$label,pred_glm)
#add model rf prediction score to class_result and ae_predict
a_glm<-as.data.frame(t(c(glm_accuracy,glm_AUC,glm_F1,glm_PRAUC,glm_MCC,glm_sensitivity,glm_specific,glm_classloss)))
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
train_label<-RFtrain$label
combo_train<-data.frame(svm_pred_train,rf_pred_train,nn_pred_train,lvq_pred_train,gbm_pred_train,glm_pred_train,train_label)
levels(combo_train$train_label)<-make.names(levels(factor(combo_train$train_label)))
label<-as.factor(test_MM$label)
combo_test<-data.frame(pred,pred_rf,pred_nn,pred_lvq,pred_gbm,pred_glm,label)
colnames(combo_test)<-colnames(combo_train) 
levels(combo_test$train_label)<-make.names(levels(factor(combo_test$train_label)))
fit_stacked<-caret::train(as.factor(train_label)~.,
                          combo_train,
                          method="nnet",
                          tuneGrid=nn_Grid,
                          MaxNWts=10000,
                          maxit=100,
                          trace=FALSE,
                          metric="Accuracy",
                          trControl=stacked_tc)
prob_stacked<-predict(fit_stacked,combo_test[,1:7],type='prob')
th<-0.5
pred_stacked<-factor(ifelse(prob_stacked$X1>th,"1","-1"))
pred_stacked<-as.vector(pred_stacked)
test_MM$label<-as.vector(test_MM$label)
stacked_accuracy<-Accuracy(pred_stacked,test_MM$label)
stacked_AUC<-adjustAUC(pred_stacked,as.numeric(test_MM$label))
stacked_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_stacked))
stacked_PRAUC<-PRAUC(as.numeric(pred_stacked),as.numeric(test_MM$label))
stacked_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_stacked),cutoff=1)
stacked_sensitivity<-Sensitivity(test_MM$label,pred_stacked)
stacked_specific<-Specificity(test_MM$label,pred_stacked)
stacked_classloss<-ZeroOneLoss(test_MM$label,pred_stacked)
#add model rf prediction score to class_result and ae_predict
a_stacked<-as.data.frame(t(c(stacked_accuracy,stacked_AUC,stacked_F1,stacked_PRAUC,stacked_MCC,stacked_sensitivity,stacked_specific,stacked_classloss)))
a_stacked<-as.data.frame(round(a_stacked,4))
a_stacked.model<-paste("stacked_model",sep="")
a_stacked<-cbind(a_stacked.model,a_stacked)
colnames(a_stacked)<-colnames(ae_predict)
ae_predict<-rbind(ae_predict,a_stacked)
pred_stacked.model<-paste("stacked_model",sep="")
pred_stacked_class<-cbind(pred_stacked.model,t(pred_stacked))
colnames(pred_stacked_class)<-colnames(class_result)
class_result<-rbind(class_result,pred_stacked_class)

#save model seperately
saveRDS(svmfit,"./II_svmfit.rds")
saveRDS(rffit,"./II_rffit.rds")
saveRDS(nn_fit,"./II_nn_fit.rds")
saveRDS(lvq_fit,"./II_lvq_fit.rds")
saveRDS(gbm_fit,"./II_gbm_fit.rds")
saveRDS(glmfit,"./II_glmfit.rds")
saveRDS(fit_stacked,"./II_fit_stacked.rds")



####ISS=III
train_ISS_III<-subset(train_data,D_ISS=="3")
b<-c("D_ISS")
train_ISS_III<-train_ISS_III[,!(names(train_ISS_III)%in%b)]
data_III<-ubBalance(X=train_ISS_III[,1:(ncol(train_ISS_III)-1)], Y=train_ISS_III$label, type="ubSMOTE", positive=1, percOver=100, percUnder=200,
                  k=5, perc=50, method="percPos", w=NULL, verbose=FALSE)
data_III<-cbind(data_III$X,data_III$Y)
train_III<-data_III
colnames(train_III)[ncol(train_III)]<-c("label")
set.seed(123)
train.rows<-createDataPartition(y=train_III$label,p=0.85,list=FALSE)
train_MM<-train_III[train.rows,]
test_MM<-train_III[-train.rows,]

nfea<-ncol(train_MM)-1




##create a dataframe to store prediction score for different model
ae_predict<-setNames(as.data.frame(matrix(ncol=9,nrow = 0)),c("model","Accuracy","AUC","F1_score","PRAUC","MCC","sensitivity","specific","classloss"))
class_result<-as.data.frame(matrix(ncol=nrow(test_MM)+1,nrow = 0))
pred_score<-as.data.frame(matrix(ncol=nrow(test_MM)+1,nrow = 0))
#function change label to 0 and 1 for calculate AUC and MMC
adjustAUC<-function(predict,truth){
  {for (i in 1:length(predict)) {if (predict[i]==-1 ){{predict[i]=0}}}}
  {for (i in 1:length(truth)) {if (truth[i]==-1 ){{truth[i]=0}}}}
  return(AUC(predict,truth))
}

adjustmcc<-function(truth,predict,cutoff=1){
  {for (i in 1:length(predict)) {if (predict[i]==-1 ){{predict[i]=0}}}}
  {for (i in 1:length(truth)) {if (truth[i]==-1 ){{truth[i]=0}}}}
  return(mcc(truth,predict,cutoff=1))
}


##SVM
##the random number seed is set before each algorithm is trained to ensure that each algorithm gets the same data partitions and repeats
library(e1071)
tc<-tune.control(cross=10)
svmtrain<-train_MM
svmfit<-svm(x=svmtrain[,1:nfea],y=as.factor(svmtrain$label),gamma=0.01,tune.control=tc,scale=TRUE,kernel = 'radial',probability=TRUE)
svm_pred_train<-predict(svmfit,svmtrain[,1:nfea])
prob<-predict(svmfit,test_MM[,1:nfea],probability=TRUE)
pred<-predict(svmfit,test_MM[,1:nfea])
pred<-as.vector(pred)
test_MM$label<-as.vector(test_MM$label)
SVM_accuracy<-Accuracy(pred,test_MM$label)
SVM_AUC<-adjustAUC(as.numeric(pred),as.numeric(test_MM$label))
SVM_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred))
SVM_PRAUC<-PRAUC(as.numeric(pred),as.numeric(test_MM$label))
SVM_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred),cutoff=1)
SVM_sensitivity<-Sensitivity(test_MM$label,pred)
SVM_specific<-Specificity(test_MM$label,pred)
SVM_classloss<-ZeroOneLoss(test_MM$label,pred)
#add model svm prediction score to class_result and ae_predict
asvm<-as.data.frame(t(c(SVM_accuracy,SVM_AUC,SVM_F1,SVM_PRAUC,SVM_MCC,SVM_sensitivity,SVM_specific,SVM_classloss)))
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
levels(RFtrain$label)<-make.names(levels(factor(RFtrain$label)))
rffit<-caret::train(x=RFtrain[,1:nfea],
                    y=RFtrain$label,
                    method="parRF",
                    tuneGrid=expand.grid(
                      .mtry=c(5,3,2,1)),
                    metric="Accuracy",
                    preProcess = c('center', 'scale'),
                    exclude=NULL,
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
rf_AUC<-adjustAUC(as.numeric(pred_rf),as.numeric(test_MM$label))
rf_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_rf))
rf_PRAUC<-PRAUC(as.numeric(pred_rf),as.numeric(test_MM$label))
rf_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_rf),cutoff=1)
rf_sensitivity<-Sensitivity(test_MM$label,pred_rf)
rf_specific<-Specificity(test_MM$label,pred_rf)
rf_classloss<-ZeroOneLoss(test_MM$label,pred_rf)
#add model rf prediction score to class_result and ae_predict
arf<-as.data.frame(t(c(rf_accuracy,rf_AUC,rf_F1,rf_PRAUC,rf_MCC,rf_sensitivity,rf_specific,rf_classloss)))
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
  .size=c(30,5),
  .decay=0.00147)
nn_fit<-caret::train(label~.,
                     NNtrain,
                     method="nnet",
                     metric="Accuracy",
                     preProcess = c('center', 'scale'),
                     tuneGrid=nn_Grid,
                     MaxNWts=10000,
                     maxit=50,
                     trace=FALSE)
pred_nn<-predict(nn_fit,test_MM[,1:nfea])
nn_pred_train<-predict(nn_fit,RFtrain[,1:nfea])
prob_nn<-predict(nn_fit,test_MM,type='prob')
nn_prob_train<-predict(nn_fit,RFtrain[,1:nfea],type="prob")
nn_pred_train<-predict(nn_fit,RFtrain[,1:nfea])
#prob_nn<-predict(nn_fit,test_MM,type='prob')
pred_nn<-as.vector(pred_nn)
test_MM$label<-as.vector(test_MM$label)
nn_accuracy<-Accuracy(pred_nn,test_MM$label)
nn_AUC<-adjustAUC(as.numeric(pred_nn),as.numeric(test_MM$label))
nn_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_nn))
nn_PRAUC<-PRAUC(as.numeric(pred_nn),as.numeric(test_MM$label))
nn_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_nn),cutoff=1)
nn_sensitivity<-Sensitivity(test_MM$label,pred_nn)
nn_specific<-Specificity(test_MM$label,pred_nn)
nn_classloss<-ZeroOneLoss(test_MM$label,pred_nn)
#add model rf prediction score to class_result and ae_predict
a_nn<-as.data.frame(t(c(nn_accuracy,nn_AUC,nn_F1,nn_PRAUC,nn_MCC,nn_sensitivity,nn_specific,nn_classloss)))
a_nn<-as.data.frame(round(a_nn,4))
a_nn.model<-paste("ANN",sep="")
a_nn<-cbind(a_nn.model,a_nn)
colnames(a_nn)<-colnames(ae_predict)
ae_predict<-rbind(ae_predict,a_nn)
pred_nn.model<-paste("ANN",sep="")
pred_nn_class<-cbind(pred_nn.model,t(pred_nn))
colnames(pred_nn_class)<-colnames(class_result)
class_result<-rbind(class_result,pred_nn_class)
#Learning Vector Quantization -> xgboost
set.seed(7)
lvqcontrol<-trainControl(method = "cv", number = 3, returnResamp = "all", search = "random",classProbs = TRUE, summaryFunction = twoClassSummary)
LVQtrain<-train_MM
lvq_grid <-expand.grid(nrounds = c(1, 10),
                       max_depth = c(1, 4),
                       eta = c(.1, .4),
                       gamma = 0,
                       colsample_bytree = .7,
                       min_child_weight = 1,
                       subsample = c(.8, 1))
lvq_fit<-caret::train(label~.,
                      RFtrain,
                      method="xgbTree",
                      trControl=lvqcontrol,
                      preProcess = c('center', 'scale'),
                      turnGrid=lvq_grid,
                      metric="ROC",
                      verbose=FALSE)
prob_lvq<-predict(lvq_fit,test_MM[,1:nfea],type='prob')
th<-0.5
pred_lvq<-factor(ifelse(prob_lvq$X1>th,"1","-1"))
lvq_prob_train<-predict(lvq_fit,RFtrain[,1:nfea],type="prob")
lvq_pred_train<-factor(ifelse(lvq_prob_train$X1>th,"1","-1"))
pred_lvq<-as.vector(pred_lvq)
test_MM$label<-as.vector(test_MM$label)
lvq_accuracy<-Accuracy(pred_lvq,test_MM$label)
lvq_AUC<-adjustAUC(pred_lvq,as.numeric(test_MM$label))
lvq_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_lvq))
lvq_PRAUC<-PRAUC(as.numeric(pred_lvq),as.numeric(test_MM$label))
lvq_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_lvq),cutoff=1)
lvq_sensitivity<-Sensitivity(test_MM$label,pred_lvq)
lvq_specific<-Specificity(test_MM$label,pred_lvq)
lvq_classloss<-ZeroOneLoss(test_MM$label,pred_lvq)
#add model lvq prediction score to class_result and ae_predict
a_lvq<-as.data.frame(t(c(lvq_accuracy,lvq_AUC,lvq_F1,lvq_PRAUC,lvq_MCC,lvq_sensitivity,lvq_specific,lvq_classloss)))
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
gbm_fit<-caret::train(label~.,
                      RFtrain,
                      method="gbm",
                      preProcess = c('center', 'scale'),
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
gbm_AUC<-adjustAUC(pred_gbm,as.numeric(test_MM$label))
gbm_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_gbm))
gbm_PRAUC<-PRAUC(as.numeric(pred_gbm),as.numeric(test_MM$label))
gbm_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_gbm),cutoff=1)
gbm_sensitivity<-Sensitivity(test_MM$label,pred_gbm)
gbm_specific<-Specificity(test_MM$label,pred_gbm)
gbm_classloss<-ZeroOneLoss(test_MM$label,pred_gbm)
#add model gbm prediction score to class_result and ae_predict
a_gbm<-as.data.frame(t(c(gbm_accuracy,gbm_AUC,gbm_F1,gbm_PRAUC,gbm_MCC,gbm_sensitivity,gbm_specific,gbm_classloss)))
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
gridsearch_for_lambda =  data.frame (alpha =  seq(.05, 1, length = 15),
                                     lambda = c ((1:5)/10))
train_control = trainControl (method = "LOOCV", returnResamp = "all", search = "random",classProbs = TRUE,summaryFunction = twoClassSummary )
glmfit<-caret::train(label~.,
                     RFtrain,
                     method="glmnet",
                     tuneGrid=gridsearch_for_lambda,
                     preProcess = c('center', 'scale'),
                     metric="ROC",
                     trControl=train_control)
prob_glm<-predict(glmfit,test_MM[,1:nfea],type='prob')
th<-0.5
pred_glm<-factor(ifelse(prob_glm$X1>th,"1","-1"))
glm_prob_train<-predict(glmfit,RFtrain[,1:nfea],type="prob")
glm_pred_train<-factor(ifelse(glm_prob_train$X1>th,"1","-1"))
pred_glm<-as.vector(pred_glm)
test_MM$label<-as.vector(test_MM$label)
glm_accuracy<-Accuracy(pred_glm,test_MM$label)
glm_AUC<-adjustAUC(pred_glm,as.numeric(test_MM$label))
glm_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_glm))
glm_PRAUC<-PRAUC(as.numeric(pred_glm),as.numeric(test_MM$label))
glm_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_glm),cutoff=1)
glm_sensitivity<-Sensitivity(test_MM$label,pred_glm)
glm_specific<-Specificity(test_MM$label,pred_glm)
glm_classloss<-ZeroOneLoss(test_MM$label,pred_glm)
#add model rf prediction score to class_result and ae_predict
a_glm<-as.data.frame(t(c(glm_accuracy,glm_AUC,glm_F1,glm_PRAUC,glm_MCC,glm_sensitivity,glm_specific,glm_classloss)))
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
train_label<-RFtrain$label
combo_train<-data.frame(svm_pred_train,rf_pred_train,nn_pred_train,lvq_pred_train,gbm_pred_train,glm_pred_train,train_label)
levels(combo_train$train_label)<-make.names(levels(factor(combo_train$train_label)))
label<-as.factor(test_MM$label)
combo_test<-data.frame(pred,pred_rf,pred_nn,pred_lvq,pred_gbm,pred_glm,label)
colnames(combo_test)<-colnames(combo_train) 
levels(combo_test$train_label)<-make.names(levels(factor(combo_test$train_label)))
fit_stacked<-caret::train(as.factor(train_label)~.,
                          combo_train,
                          method="parRF",
                          tuneGrid=expand.grid(
                            .mtry=c(3,2,1)),
                          metric="Accuracy",
                          trControl=stacked_tc)
prob_stacked<-predict(fit_stacked,combo_test[,1:7],type='prob')
th<-0.5
pred_stacked<-factor(ifelse(prob_stacked$X1>th,"1","-1"))
pred_stacked<-as.vector(pred_stacked)
test_MM$label<-as.vector(test_MM$label)
stacked_accuracy<-Accuracy(pred_stacked,test_MM$label)
stacked_AUC<-adjustAUC(pred_stacked,as.numeric(test_MM$label))
stacked_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_stacked))
stacked_PRAUC<-PRAUC(as.numeric(pred_stacked),as.numeric(test_MM$label))
stacked_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_stacked),cutoff=1)
stacked_sensitivity<-Sensitivity(test_MM$label,pred_stacked)
stacked_specific<-Specificity(test_MM$label,pred_stacked)
stacked_classloss<-ZeroOneLoss(test_MM$label,pred_stacked)
#add model rf prediction score to class_result and ae_predict
a_stacked<-as.data.frame(t(c(stacked_accuracy,stacked_AUC,stacked_F1,stacked_PRAUC,stacked_MCC,stacked_sensitivity,stacked_specific,stacked_classloss)))
a_stacked<-as.data.frame(round(a_stacked,4))
a_stacked.model<-paste("stacked_model",sep="")
a_stacked<-cbind(a_stacked.model,a_stacked)
colnames(a_stacked)<-colnames(ae_predict)
ae_predict<-rbind(ae_predict,a_stacked)
pred_stacked.model<-paste("stacked_model",sep="")
pred_stacked_class<-cbind(pred_stacked.model,t(pred_stacked))
colnames(pred_stacked_class)<-colnames(class_result)
class_result<-rbind(class_result,pred_stacked_class)

#save model seperately
saveRDS(svmfit,"./III_svmfit.rds")
saveRDS(rffit,"./III_rffit.rds")
saveRDS(nn_fit,"./III_nn_fit.rds")
saveRDS(lvq_fit,"./III_lvq_fit.rds")
saveRDS(gbm_fit,"./III_gbm_fit.rds")
saveRDS(glmfit,"./III_glmfit.rds")
saveRDS(fit_stacked,"./III_fit_stacked.rds")



I_svmfit<-readRDS("./I_svmfit.rds")
I_rffit<-readRDS("./I_rffit.rds")
I_nn_fit<-readRDS("./I_nn_fit.rds")
I_lvq_fit<-readRDS("./I_lvq_fit.rds")
I_gbm_fit<-readRDS("./I_gbm_fit.rds")
I_glmfit<-readRDS("./I_glmfit.rds")
I_fit_stacked<-readRDS("./I_fit_stacked.rds")



II_svmfit<-readRDS("./II_svmfit.rds")
II_rffit<-readRDS("./II_rffit.rds")
II_nn_fit<-readRDS("./II_nn_fit.rds")
II_lvq_fit<-readRDS("./II_lvq_fit.rds")
II_gbm_fit<-readRDS("./II_gbm_fit.rds")
II_glmfit<-readRDS("./II_glmfit.rds")
II_fit_stacked<-readRDS("./II_fit_stacked.rds")



III_svmfit<-readRDS("./III_svmfit.rds")
III_rffit<-readRDS("./III_rffit.rds")
III_nn_fit<-readRDS("./III_nn_fit.rds")
III_lvq_fit<-readRDS("./III_lvq_fit.rds")
III_gbm_fit<-readRDS("./III_gbm_fit.rds")
III_glmfit<-readRDS("./III_glmfit.rds")
III_fit_stacked<-readRDS("./III_fit_stacked.rds")

M_genes_eu_list<-model.entrezs

save(I_svmfit,I_rffit,I_nn_fit,I_lvq_fit,I_gbm_fit,I_glmfit,I_fit_stacked,
     II_svmfit,II_rffit,II_nn_fit,II_lvq_fit,II_gbm_fit,II_glmfit,II_fit_stacked,
     III_svmfit,III_rffit,III_nn_fit,III_lvq_fit,III_gbm_fit,III_glmfit,III_fit_stacked,
     M_genes_eu_list,file="model-state-metadata_ISS.Rd")
