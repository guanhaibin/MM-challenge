install.packages("MESS") #auc
library(MESS)


#Read data in R and preparation
setwd("U:/hguan003/MM")
Gene_expression=read.csv("GSE24080UAMSentrezIDlevel.csv",header=T,sep=",",row.names=1)
clinical<-read.csv("globalClinTraining.csv",header=T,sep=",",row.names = 2,stringsAsFactors = FALSE)
a<-c("GSE24080UAMS","EMTAB4032","HOVON65")
clinical<-subset(clinical,Study%in%a)
clinical<-subset(clinical[,c("Study","D_Age","D_OS","D_OS_FLAG","D_PFS","D_PFS_FLAG","D_ISS","HR_FLAG")])
clinical$Patient<-rownames(clinical)

#Generate risk score RS based on D_PFS and D_PFS_FLAG
a<-subset(clinical,D_PFS_FLAG=="0")
b<-subset(clinical,D_PFS_FLAG=="1")
a_range<-as.numeric(summary(a$D_PFS)[6]-summary(a$D_PFS)[1])
b_range<-as.numeric(summary(b$D_PFS)[6]-summary(b$D_PFS)[1])
k<-5/as.numeric(-a_range)
a$RS=k*a$D_PFS+5
b$RS=k*b$D_PFS+5-k*(b_range)
RS_18=k*18*30.5+5-k*(b_range) #RS_18=9.07777
clinical_RS<-rbind(a,b)
clinical$RS<-clinical_RS$RS[match(rownames(clinical),rownames(clinical_RS))]
clinical$D_Age<-clinical$D_Age/100
clinical$D_ISS<-clinical$D_ISS/10
clinical<-subset(clinical[,c("Patient","Study","D_Age","D_ISS","RS","HR_FLAG")])
clinical_UAMS<-subset(clinical,Study=="GSE24080UAMS")

train_data=as.data.frame(t(Gene_expression))
train_data<-as.data.frame(scale(train_data))
train_data<-merge(train_data,clinical_UAMS[,c("D_Age","D_ISS","RS","HR_FLAG")],by="row.names",all.x=TRUE)
row.names(train_data)<-train_data[,1]
colnames(train_data)[1]<-c("Patient")



##Split to training/testing set
library(h2o)
h2o.init(nthreads=-1,max_mem_size = "16G",enable_assertions = FALSE)
total.hex<-as.h2o(train_data[,2:(ncol(train_data)-4)])
total.hex$label<-as.h2o(train_data[,(ncol(train_data)-1)])
splits<-h2o.splitFrame(total.hex,c(0.6,0.2),destination_frames=c("train","valid","test"),seed=1234)
train<-h2o.assign(splits[[1]],"train.hex")#60%
valid<-h2o.assign(splits[[2]],"valid.hex")#20%
test<-h2o.assign(splits[[3]],"test.hex")#20%

total_label<-as.vector(total.hex[,ncol(total.hex)])
train_label<-as.vector(train[,ncol(train)])
train[,ncol(train)]<-NULL
valid_label<-as.vector(valid[,ncol(valid)])
valid[,ncol(valid)]<-NULL
test_label<-as.vector(test[,ncol(test)])
test[,ncol(test)]<-NULL
total.hex[,ncol(total.hex)]<-NULL



hyper_params_ae<-list(
  hidden=list(c(1000,100,1000),c(500)))


ae_grid<-h2o.grid(
  algorithm="deeplearning",
  grid_id = "ae_grid_id",
  training_frame=total.hex,
  epochs=10,
  export_weights_and_biases = T,
  ignore_const_cols = F,
  autoencoder = T,
  activation=c("Tanh"),
  l1=1e-5,
  l2=1e-5,
  max_w2=10,
  variable_importances = T,
  hyper_params = hyper_params_ae
)

summary(ae_grid)
nmodel<-nrow(ae_grid@summary_table)

##create a dataframe to store prediction score for different model
ae_predict<-setNames(as.data.frame(matrix(ncol=9,nrow = 0)),c("model","AUC","PRAUC","Accuracy","F1_score","MCC","sensitivity","specific","classloss"))
Risk_Score<-as.data.frame(matrix(ncol=nrow(test)+1,nrow = 0))
class_result<-as.data.frame(matrix(ncol=nrow(test)+1,nrow = 0))

#function change label to 0 and 1 for calculate MMC
adjustmcc<-function(truth,predict,cutoff=1){
  {for (i in 1:length(predict)) {if (predict[i]==-1 ){{predict[i]=0}}}}
  {for (i in 1:length(truth)) {if (truth[i]==-1 ){{truth[i]=0}}}}
  return(mcc(truth,predict,cutoff=1))
}


for (i in 1:nmodel) {
  model<-h2o.getModel(ae_grid@model_ids[[i]])
  fealayer<-(length(model@parameters$hidden)+1)/2
  nfea<-model@parameters$hidden[fealayer]
  deep.fea<-as.data.frame(h2o.deepfeatures( model,total.hex,layer=fealayer))
  deep.fea$label<-as.character(total_label)
  deep.fea.train<-as.data.frame(h2o.deepfeatures( model,train,layer=fealayer))
  deep.fea.train$label<-as.character(train_label)
  deep.fea.valid<-as.data.frame(h2o.deepfeatures(model,valid,layer=fealayer))
  deep.fea.valid$label<-as.character(valid_label)
  deep.fea.test<-as.data.frame(h2o.deepfeatures(model,test,layer=fealayer))
  deep.fea.test$label<-as.character(test_label)
  ##SVM
  library(e1071)
  tc<-tune.control(cross=5)
  svmtrain<-rbind(deep.fea.train,deep.fea.valid)
  svmfit<-svm(svmtrain[,1:nfea],as.factor(svmtrain$label),tune.control=tc)
  pred<-predict(svmfit,deep.fea.test[,1:nfea])
  pred<-as.vector(pred)
  deep.fea.test$label<-as.vector(deep.fea.test$label)
  SVM_AUC<-MLmetrics::AUC(as.numeric(pred),as.numeric(deep.fea.test$label))
  SVM_PRAUC<-PRAUC(as.numeric(pred),as.numeric(deep.fea.test$label))
  SVM_accuracy<-Accuracy(pred,deep.fea.test$label)
  SVM_F1<-F1_Score(as.numeric(deep.fea.test$label),as.numeric(pred))
  SVM_MCC<-adjustmcc(as.numeric(deep.fea.test$label),as.numeric(pred),cutoff=1)
  SVM_sensitivity<-Sensitivity(deep.fea.test$label,pred)
  SVM_specific<-Specificity(deep.fea.test$label,pred)
  SVM_classloss<-ZeroOneLoss(deep.fea.test$label,pred)
  #add model svm prediction score to class_result and ae_predict
  asvm<-as.data.frame(t(c(SVM_AUC,SVM_PRAUC,SVM_accuracy,SVM_F1,SVM_MCC,SVM_sensitivity,SVM_specific,SVM_classloss)))
  asvm<-as.data.frame(round(asvm,4))
  asvm.model<-paste(ae_grid@summary_table[i,1],"SVM",sep="")
  asvm<-cbind(asvm.model,asvm)
  colnames(asvm)<-colnames(ae_predict)
  ae_predict<-rbind(ae_predict,asvm)
  pred.model<-paste(ae_grid@summary_table[i,1],"SVM",sep="")
  pred_class<-cbind(pred.model,t(pred))
  colnames(pred_class)<-colnames(class_result)
  class_result<-rbind(class_result,pred_class)
  ##random forest
  library(ranger)
  fitcontrol<-trainControl(method="repeatedcv",number=10,repeats=2)
  RFtrain<-rbind(deep.fea.train,deep.fea.valid)
  rffit<-caret:::train(label~.,
                       RFtrain,
                       method="ranger",
                       tuneGrid=expand.grid(
                         .mtry=2),
                       metric="Accuracy",
                       trControl=fitcontrol)
  pred_rf<-predict(rffit,deep.fea.test)
  pred_rf<-as.vector(pred_rf)
  deep.fea.test$label<-as.vector(deep.fea.test$label)
  rf_AUC<-AUC(as.numeric(pred_rf),as.numeric(deep.fea.test$label))
  rf_accuracy<-Accuracy(pred_rf,deep.fea.test$label)
  rf_F1<-F1_Score(as.numeric(deep.fea.test$label),as.numeric(pred_rf))
  rf_PRAUC<-PRAUC(as.numeric(pred_rf),as.numeric(deep.fea.test$label))
  rf_MCC<-adjustmcc(as.numeric(deep.fea.test$label),as.numeric(pred_rf),cutoff=1)
  rf_sensitivity<-Sensitivity(deep.fea.test$label,pred_rf)
  rf_specific<-Specificity(deep.fea.test$label,pred_rf)
  rf_classloss<-ZeroOneLoss(deep.fea.test$label,pred_rf)
  #add model rf prediction score to class_result and ae_predict
  arf<-as.data.frame(t(c(rf_AUC,rf_PRAUC,rf_accuracy,rf_F1,rf_MCC,rf_sensitivity,rf_specific,rf_classloss)))
  arf<-as.data.frame(round(arf,4))
  arf.model<-paste(ae_grid@summary_table[i,1],"RF",sep="")
  arf<-cbind(arf.model,arf)
  colnames(arf)<-colnames(ae_predict)
  ae_predict<-rbind(ae_predict,arf)
  pred_rf.model<-paste(ae_grid@summary_table[i,1],"RF",sep="")
  pred_rf_class<-cbind(pred_rf.model,t(pred_rf))
  colnames(pred_rf_class)<-colnames(class_result)
  class_result<-rbind(class_result,pred_rf_class)
  ##K nearest neighbor
  knnfit<-caret::train(label~.,
                       RFtrain,
                       method="knn",
                       tuneGrid=expand.grid(
                         .k=5),
                       trControl=fitcontrol)
  pred_knn<-predict(knnfit,deep.fea.test)
  pred_knn<-as.vector(pred_knn)
  deep.fea.test$label<-as.vector(deep.fea.test$label)
  knn_accuracy<-Accuracy(pred_knn,deep.fea.test$label)
  knn_AUC<-AUC(as.numeric(pred_knn),as.numeric(deep.fea.test$label))
  knn_F1<-F1_Score(as.numeric(deep.fea.test$label),as.numeric(pred_knn))
  knn_PRAUC<-PRAUC(as.numeric(pred_knn),as.numeric(deep.fea.test$label))
  knn_MCC<-adjustmcc(as.numeric(deep.fea.test$label),as.numeric(pred_knn),cutoff=1)
  knn_sensitivity<-Sensitivity(deep.fea.test$label,pred_knn)
  knn_specific<-Specificity(deep.fea.test$label,pred_knn)
  knn_classloss<-ZeroOneLoss(deep.fea.test$label,pred_knn)
  #add model rf prediction score to class_result and ae_predict
  aknn<-as.data.frame(t(c(knn_AUC,knn_PRAUC,knn_accuracy,knn_F1,knn_MCC,knn_sensitivity,knn_specific,knn_classloss)))
  aknn<-as.data.frame(round(aknn,4))
  aknn.model<-paste(ae_grid@summary_table[i,1],"KNN",sep="")
  aknn<-cbind(aknn.model,aknn)
  colnames(aknn)<-colnames(ae_predict)
  ae_predict<-rbind(ae_predict,aknn)
  pred_knn.model<-paste(ae_grid@summary_table[i,1],"KNN",sep="")
  pred_knn_class<-cbind(pred_knn.model,t(pred_knn))
  colnames(pred_knn_class)<-colnames(class_result)
  class_result<-rbind(class_result,pred_knn_class)
  #Multiple Paerceptron Network by Stachastic Gradient Descent
  nn_Grid<-expand.grid(
    .size=c(50,10,5),
    .decay=0.00147)
  nn_fit<-caret::train(label~.,
                       RFtrain,
                       method="nnet",
                       metric="Accuracy",
                       tuneGrid=nn_Grid,
                       MaxNWts=10000,
                       maxit=100,
                       trControl=fitcontrol,
                       trace=FALSE)
  pred_nn<-predict(nn_fit,deep.fea.test)
  pred_nn<-as.vector(pred_nn)
  deep.fea.test$label<-as.vector(deep.fea.test$label)
  nn_accuracy<-Accuracy(pred_nn,deep.fea.test$label)
  nn_AUC<-AUC(as.numeric(pred_nn),as.numeric(deep.fea.test$label))
  nn_F1<-F1_Score(as.numeric(deep.fea.test$label),as.numeric(pred_nn))
  nn_PRAUC<-PRAUC(as.numeric(pred_nn),as.numeric(deep.fea.test$label))
  nn_MCC<-adjustmcc(as.numeric(deep.fea.test$label),as.numeric(pred_nn),cutoff=1)
  nn_sensitivity<-Sensitivity(deep.fea.test$label,pred_nn)
  nn_specific<-Specificity(deep.fea.test$label,pred_nn)
  nn_classloss<-ZeroOneLoss(deep.fea.test$label,pred_nn)
  #add model rf prediction score to class_result and ae_predict
  a_nn<-as.data.frame(t(c(nn_AUC,nn_PRAUC,nn_accuracy,nn_F1,nn_MCC,nn_sensitivity,nn_specific,nn_classloss)))
  a_nn<-as.data.frame(round(a_nn,4))
  a_nn.model<-paste(ae_grid@summary_table[i,1],"ANN",sep="")
  a_nn<-cbind(a_nn.model,a_nn)
  colnames(a_nn)<-colnames(ae_predict)
  ae_predict<-rbind(ae_predict,a_nn)
  pred_nn.model<-paste(ae_grid@summary_table[i,1],"ANN",sep="")
  pred_nn_class<-cbind(pred_nn.model,t(pred_nn))
  colnames(pred_nn_class)<-colnames(class_result)
  class_result<-rbind(class_result,pred_nn_class)
  ##Stacked model
  #combine all the predictions of above classifiers
  combo<-data.frame(pred,pred_rf,pred_knn,pred_nn,label=deep.fea.test$label)
  fit_stacked<-caret::train(as.factor(label)~.,
                            combo,
                            method="ranger",
                            metric="Accuracy")
  pred_stacked<-predict(fit_stacked,deep.fea.test$label)
  pred_stacked<-as.vector(pred_stacked)
  deep.fea.test$label<-as.vector(deep.fea.test$label)
  stacked_accuracy<-Accuracy(pred_stacked,deep.fea.test$label)
  stacked_AUC<-AUC(as.numeric(pred_stacked),as.numeric(deep.fea.test$label))
  stacked_F1<-F1_Score(as.numeric(deep.fea.test$label),as.numeric(pred_stacked))
  stacked_PRAUC<-PRAUC(as.numeric(pred_stacked),as.numeric(deep.fea.test$label))
  stacked_MCC<-adjustmcc(as.numeric(deep.fea.test$label),as.numeric(pred_stacked),cutoff=1)
  stacked_sensitivity<-Sensitivity(deep.fea.test$label,pred_stacked)
  stacked_specific<-Specificity(deep.fea.test$label,pred_stacked)
  stacked_classloss<-ZeroOneLoss(deep.fea.test$label,pred_stacked)
  #add model rf prediction score to class_result and ae_predict
  a_stacked<-as.data.frame(t(c(stacked_AUC,stacked_PRAUC,stacked_accuracy,stacked_F1,stacked_MCC,stacked_sensitivity,stacked_specific,stacked_classloss)))
  a_stacked<-as.data.frame(round(a_stacked,4))
  a_stacked.model<-paste(ae_grid@summary_table[i,1],"stacked_model",sep="")
  a_stacked<-cbind(a_stacked.model,a_stacked)
  colnames(a_stacked)<-colnames(ae_predict)
  ae_predict<-rbind(ae_predict,a_stacked)
  pred_stacked.model<-paste(ae_grid@summary_table[i,1],"stacked_model",sep="")
  pred_stacked<-cbind(pred_stacked.model,t(pred_stacked))
  colnames(pred_stacked)<-colnames(class_result)
  class_result<-rbind(class_result,pred_stacked)
}


write.csv(class_result,file="HR_FLAG.csv")
write.csv(ae_predict,file="model_predict_score.csv")