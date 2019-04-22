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


#read UAMS datasets with selected genes
Gene_expression=read.csv("GSE24080UAMSentrezIDlevel.csv",header=T,sep=",",row.names=1)
clinical<-read.csv("globalClinTraining.csv",header=T,sep=",",row.names = 2,stringsAsFactors = FALSE)
clinical<-subset(clinical,Study=="GSE24080UAMS")
clinical<-subset(clinical[,c("Study","HR_FLAG")])
clinical[clinical[,"HR_FLAG"]==TRUE,2]=1
clinical[clinical[,"HR_FLAG"]==FALSE,2]=-1
train_data<-as.data.frame(t(scale(Gene_expression)))
train_data<-merge(train_data,clinical,by="row.names",all.x=TRUE)
row.names(train_data)<-train_data[,1]
train_data<-train_data[-c(1)]
colnames(train_data)[ncol(train_data)]<-c("label")
UAMS<-train_data
UAMS<-UAMS[UAMS$label!="CENSORED",]
UAMS$label<-as.factor(UAMS$label)
UAMS_test<-UAMS[,!(colnames(UAMS)%in%c("Study"))]

#split UAMS test data into different subset to fit in different models
UAMS_U_test<-subset(UAMS_test,select=c(U_genes_list))
UAMS_U_test$label<-UAMS_test$label
UAMS_E_test<-setNames(as.data.frame(matrix(0,nrow=nrow(UAMS_test),ncol=length(E_genes_list))),E_genes_list)
idx<-colnames(UAMS_test)[colnames(UAMS_test)%in%colnames(UAMS_E_test)]
UAMS_E_test[,idx]<-UAMS_test[,idx]
rownames(UAMS_E_test)<-rownames(UAMS_test)
UAMS_E_test$label<-UAMS_test$label
UAMS_H_test<-subset(UAMS_test,select=c(H_genes_list))
UAMS_H_test$label<-UAMS_test$label

#read EMTAB datasets with selected genes
Gene_expression=read.csv("EMTAB4032entrezIDlevel.csv",header=T,sep=",",row.names=1)
clinical<-read.csv("globalClinTraining.csv",header=T,sep=",",row.names = 2,stringsAsFactors = FALSE)
clinical<-subset(clinical,Study=="EMTAB4032")
colnames(Gene_expression)<-str_replace_all(colnames(Gene_expression),"\\.","-")
clinical<-clinical[setdiff(rownames(clinical),c("MM-11-136","MM-12-161")),]
clinical<-subset(clinical[,c("Study","HR_FLAG")])
clinical[clinical[,"HR_FLAG"]==TRUE,2]=1
clinical[clinical[,"HR_FLAG"]==FALSE,2]=-1
train_data<-as.data.frame(t(scale(Gene_expression)))
train_data<-merge(train_data,clinical,by="row.names",all.x=TRUE)
row.names(train_data)<-train_data[,1]
train_data<-train_data[-c(1)]
colnames(train_data)[ncol(train_data)]<-c("label")
EMTAB<-train_data
EMTAB<-EMTAB[EMTAB$label!="CENSORED",]
EMTAB$label<-as.factor(EMTAB$label)
EMTAB_test<-EMTAB[,!(colnames(EMTAB)%in%c("Study"))]



#split EMTAB test data into different subset to fit in different models
EMTAB_E_test<-subset(EMTAB_test,select=c(E_genes_list))
EMTAB_E_test$label<-EMTAB_test$label
EMTAB_U_test<-setNames(as.data.frame(matrix(0,nrow=nrow(EMTAB_test),ncol=length(U_genes_list))),U_genes_list)
idx<-colnames(EMTAB_test)[colnames(EMTAB_test)%in%colnames(EMTAB_U_test)]
EMTAB_U_test[,idx]<-EMTAB_test[,idx]
rownames(EMTAB_U_test)<-rownames(EMTAB_test)
EMTAB_U_test$label<-EMTAB_test$label
EMTAB_H_test<-setNames(as.data.frame(matrix(0,nrow=nrow(EMTAB_test),ncol=length(H_genes_list))),H_genes_list)
idx<-colnames(EMTAB_test)[colnames(EMTAB_test)%in%colnames(EMTAB_H_test)]
EMTAB_H_test[,idx]<-EMTAB_test[,idx]
rownames(EMTAB_H_test)<-rownames(EMTAB_test)
EMTAB_H_test$label<-EMTAB_test$label

#read HOVON datasets with selected genes
Gene_expression=read.csv("GSE19784HOVON65entrezIDlevel.csv",header=T,sep=",",row.names=1)
clinical<-read.csv("globalClinTraining.csv",header=T,sep=",",row.names = 2,stringsAsFactors = FALSE)
clinical<-subset(clinical,Study=="HOVON65")
clinical<-subset(clinical[,c("Study","HR_FLAG")])
clinical[clinical[,"HR_FLAG"]==TRUE,2]=1
clinical[clinical[,"HR_FLAG"]==FALSE,2]=-1
train_data<-as.data.frame(t(scale(Gene_expression)))
train_data<-merge(train_data,clinical,by="row.names",all.x=TRUE)
row.names(train_data)<-train_data[,1]
train_data<-train_data[-c(1)]
colnames(train_data)[ncol(train_data)]<-c("label")
HOVON<-train_data
HOVON<-HOVON[HOVON$label!="CENSORED",]
HOVON$label<-as.factor(HOVON$label)
HOVON_test<-HOVON[,!(colnames(HOVON)%in%c("Study"))]

#split HOVON_test data into different subset to fit in different models
HOVON_U_test<-subset(HOVON_test,select=c(U_genes_list))
HOVON_U_test$label<-HOVON_test$label
HOVON_E_test<-setNames(as.data.frame(matrix(0,nrow=nrow(HOVON_test),ncol=length(E_genes_list))),E_genes_list)
idx<-colnames(HOVON_test)[colnames(HOVON_test)%in%colnames(HOVON_E_test)]
HOVON_E_test[,idx]<-HOVON_test[,idx]
rownames(HOVON_E_test)<-rownames(HOVON_test)
HOVON_E_test$label<-HOVON_test$label
HOVON_H_test<-subset(HOVON_test,select=c(H_genes_list))
HOVON_H_test$label<-HOVON_test$label

##combine all test data correspomding to the different model
test_U<-rbind(UAMS_U_test,EMTAB_U_test,HOVON_U_test)
test_E<-rbind(UAMS_E_test,EMTAB_E_test,HOVON_E_test)
test_H<-rbind(UAMS_H_test,EMTAB_H_test,HOVON_H_test)

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
  a<-Accuracy(ensemble_test[,i],ensemble_test$label)
  ensemble_accuracy[i,1]<-a
  rownames(ensemble_accuracy)[i]<-colnames(ensemble_test)[i]
}
ensemble_accuracy[22,1]<-ensemble_rf_accuracy
rownames(ensemble_accuracy)[22]<-"ensemble_rf"
ensemble_accuracy[23,1]<-ensemble_nn_accuracy
rownames(ensemble_accuracy)[23]<-"ensemble_nn"
rownames(ensemble_accuracy)<-gsub("train","test",rownames(ensemble_accuracy))
write.csv(ensemble_accuracy,"ensemble_accuracy.csv")



###generate prediction score
predictionscore<-array()
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
UAMS_S_P<-data.frame(UAMS$Study,rownames(UAMS))
colnames(UAMS_S_P)<-c("Study","Patient")
EMTAB_S_P<-data.frame(EMTAB$Study,rownames(EMTAB))
colnames(EMTAB_S_P)<-c("Study","Patient")
HOVON_S_P<-data.frame(HOVON$Study,rownames(HOVON))
colnames(HOVON_S_P)<-c("Study","Patient")
S_P_map<-rbind(UAMS_S_P,EMTAB_S_P,HOVON_S_P)


###Output1: max P(E1),P(E2),P(E3)
output1<-data.frame(Patient,predictionscore,highriskflag)
output<-merge(S_P_map,output1,by="Patient",all.x=TRUE)
output<-output[,c("Study","Patient","predictionscore","highriskflag")]
write.csv(output,"predictionfile.csv",row.names = FALSE)



stopImplicitCluster()



