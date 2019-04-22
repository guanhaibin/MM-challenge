
library(caret)
library(gbm)
library(mlbench)
library(MLmetrics)
library(ModelMetrics)
library(caret)
library(plyr)
library(FCNN4R)
library(matrixStats)

setwd("U:/hguan003/MM_ISS_r2")
load("emc92_uams70.Rd")
genes_list<-unique(c(emc92,uams70))
eu_genes_list<-unique(genes_list)
Gene_expression=read.csv("UAMS_training.csv",header=T,sep=",",row.names=1)
U_genes_list<-rownames(Gene_expression)[-c(1,2)]

colnames(Gene_expression)<-gsub("X","U",colnames(Gene_expression))
U_ISS<-Gene_expression[1,]
colnames(U_ISS)<-colnames(Gene_expression)
U_ISS<-as.vector(t(U_ISS))


Gene_expression=read.csv("EMTAB_training.csv",header=T,sep=",",row.names=1)
colnames(Gene_expression)<-gsub("X","E",colnames(Gene_expression))
E_ISS<-Gene_expression[1,]
colnames(E_ISS)<-colnames(Gene_expression)
E_ISS<-as.vector(t(E_ISS))
E_genes_list<-read.csv("E_genes_list.csv",row.names = 1)
E_genes_list<-as.character(E_genes_list[,1])

Gene_expression=read.csv("HOVON_training.csv",header=T,sep=",",row.names=1)
H_genes_list<-rownames(Gene_expression)[-c(1,2)]

colnames(Gene_expression)<-gsub("X","H",colnames(Gene_expression))
H_ISS<-Gene_expression[1,]
colnames(H_ISS)<-colnames(Gene_expression)
H_ISS<-as.vector(t(H_ISS))

model.entrezs<-c(eu_genes_list,U_genes_list,H_genes_list,E_genes_list)
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
UAMS_data <- as.data.frame(select.data)
UAMS_data$label<-UAMS_label
UAMS_data$ISS<-U_ISS

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
HOVON_data <- select.data
HOVON_data$label<-HOVON_label
HOVON_data$ISS<-H_ISS


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
EMTAB_data <- select.data
EMTAB_data$label<-EMTAB_label
EMTAB_data$ISS<-E_ISS


E_train_MM<-read.csv("E_train_MM.csv",row.names = 1)
E_train_MM<-as.character(E_train_MM[,1])
E_test_MM<-read.csv("E_test_MM.csv",row.names = 1)
E_test_MM<-as.character(E_test_MM[,1])
E_train_MM<-subset(EMTAB_data,rownames(EMTAB_data)%in%E_train_MM)
E_test_MM<-subset(EMTAB_data,rownames(EMTAB_data)%in%E_test_MM)



H_train_MM<-read.csv("H_train_MM.csv",row.names = 1)
H_train_MM<-as.character(H_train_MM[,1])
H_test_MM<-read.csv("H_test_MM.csv",row.names = 1)
H_test_MM<-as.character(H_test_MM[,1])
H_train_MM<-subset(HOVON_data,rownames(HOVON_data)%in%H_train_MM)
H_test_MM<-subset(HOVON_data,rownames(HOVON_data)%in%H_test_MM)


U_test_MM<-read.csv("U_test_MM.csv",row.names = 1)
U_train_MM<-read.csv("U_train_MM.csv",row.names = 1)
U_test_MM<-as.character(U_test_MM[,1])
U_train_MM<-as.character(U_train_MM[,1])
U_train_MM<-subset(UAMS_data,rownames(UAMS_data)%in%U_train_MM)
U_test_MM<-subset(UAMS_data,rownames(UAMS_data)%in%U_test_MM)

train_MM<-rbind(U_train_MM,H_train_MM,E_train_MM)
test_MM<-rbind(U_test_MM,H_test_MM,E_test_MM)

train_MM_label<-train_MM$label
train_MM$label<-NULL
test_MM_label<-test_MM$label
test_MM$label<-NULL

train_MM_ISS<-train_MM$ISS
train_MM$ISS<-NULL
test_MM_ISS<-test_MM$ISS
test_MM$ISS<-NULL

U_mean<-rowMeans(as.matrix(t(train_MM)))
U_sd<-rowSds(as.matrix(t(train_MM)))
train=sweep(train_MM,MARGIN=2,U_mean,"-")
train_MM<-sweep(train,MARGIN=2,U_sd,"/")
test=sweep(test_MM,MARGIN=2,U_mean,"-")
test_MM<-sweep(test,MARGIN=2,U_sd,"/")

test_MM$ISS<-as.numeric(test_MM_ISS)
train_MM$ISS<-as.numeric(train_MM_ISS)
train_MM$label<-train_MM_label
test_MM$label<-test_MM_label
nfea<-ncol(train_MM)-1

saveRDS(U_mean,"./U_mean.rds")
saveRDS(U_sd,"./U_sd.rds")


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
svmfit<-svm(as.factor(label)~.,data=svmtrain,gamma=0.001,tune.control=tc,kernel = 'radial',probability=TRUE)
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
rffit<-caret::train(label~.,
                    RFtrain,
                    method="parRF",
                    tuneGrid=expand.grid(
                      .mtry=c(30,10,5,2,1)),
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
                          metric="Accuracy",
                          tuneGrid=nn_Grid,
                          MaxNWts=10000,
                          maxit=100,
                          trace=FALSE)
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

write.csv(class_result,file="UAMS_TRAIN_pred_HR_FLAG.csv")
write.csv(ae_predict,file="UAMS_TRAIN_model_Accuracy.csv")

#save model seperately
saveRDS(svmfit,"./U_svmfit.rds")
saveRDS(rffit,"./U_rffit.rds")
saveRDS(nn_fit,"./U_nn_fit.rds")
saveRDS(lvq_fit,"./U_lvq_fit.rds")
saveRDS(gbm_fit,"./U_gbm_fit.rds")
saveRDS(glmfit,"./U_glmfit.rds")
saveRDS(fit_stacked,"./U_fit_stacked.rds")

U_genes_eu<-model.entrezs
write.csv(U_genes_eu,"U_genes_eu.csv")


