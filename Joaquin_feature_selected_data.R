library(limma)
library(MLmetrics)
library(ModelMetrics)
library(caret)
library(plyr)
library(FCNN4R)
library(doParallel)

Gene1=read.csv("GSE24080UAMSentrezIDlevel.csv",header=T,sep=",",row.names=1)
Gene1<-as.data.frame(t(scale(Gene1)))
Gene2=read.csv("GSE19784HOVON65entrezIDlevel.csv",header=T,sep=",",row.names=1)
Gene2<-as.data.frame(t(scale(Gene2)))
# > length(common)
# [1] 20514
Gene2=read.csv("EMTAB4032entrezIDlevel.csv",header=T,sep=",",row.names=1)
Gene2<-as.data.frame(t(scale(Gene2)))
common<-intersect(colnames(Gene1),colnames(Gene2)) 
# > length(common)
#[1] 17038




setwd("U://hguan003/MM")
Gene_expression=read.csv("HOVON_subarray.csv",header=T,sep=",",row.names=1)
Gene_expression<-as.data.frame(t(scale(Gene_expression)))
HR_FLAG<-read.csv("HOVON_target.csv",header=T,sep=",")
HR_FLAG<-t(HR_FLAG)
train_data<-merge(Gene_expression,HR_FLAG,by="row.names",all.x=TRUE)
row.names(train_data)<-train_data[,1]
train_data<-train_data[-c(1)]
colnames(train_data)[ncol(train_data)]<-c("label")
train_data$label<-as.factor(train_data$label)

set.seed(222)
random_splits<-runif(nrow(train_data))
train_MM<-train_data[random_splits<0.7,] 
test_MM<-train_data[random_splits>=0.3,]
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
library(e1071)
tc<-tune.control(cross=10)
svmtrain<-train_MM
svmfit<-svm(as.factor(label)~.,data=svmtrain,gamma=0.001,tune.control=tc,kernel = 'radial')
pred<-predict(svmfit,test_MM[,1:nfea])
#prob<-attr(pred,"probabilities")[,2]
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
library(randomForest)
fitcontrol<-trainControl(method="cv",number=10)
RFtrain<-train_MM
rffit<-caret:::train(label~.,
                     RFtrain,
                     method="parRF",
                     tuneGrid=expand.grid(
                       .mtry=30),
                     metric="Accuracy",
                     trControl=fitcontrol)
pred_rf<-predict(rffit,test_MM)
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
##K nearest neighbor
#pred_knn<-caret::knn3Train(RFtrain[,1:nfea],test_MM[,1:nfea],RFtrain[,nfea+1],k=5,prob=T)
 knnfit<-caret::train(label~.,
                       RFtrain,
                       method="knn",
                       tuneGrid=expand.grid(
                         .k=5),
                       trControl=fitcontrol)
#prob_knn<-attr(pred_knn,"probabilities")[,2]
pred_knn<-predict(knnfit,test_MM)
pred_knn<-as.vector(pred_knn)
test_MM$label<-as.vector(test_MM$label)
knn_accuracy<-Accuracy(pred_knn,test_MM$label)
knn_AUC<-adjustAUC(pred_knn,as.numeric(test_MM$label))
knn_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_knn))
knn_PRAUC<-PRAUC(as.numeric(pred_knn),as.numeric(test_MM$label))
knn_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_knn),cutoff=1)
knn_sensitivity<-Sensitivity(test_MM$label,pred_knn)
knn_specific<-Specificity(test_MM$label,pred_knn)
knn_classloss<-ZeroOneLoss(test_MM$label,pred_knn)
#add model rf prediction score to class_result and ae_predict
aknn<-as.data.frame(t(c(knn_accuracy,knn_AUC,knn_F1,knn_PRAUC,knn_MCC,knn_sensitivity,knn_specific,knn_classloss)))
aknn<-as.data.frame(round(aknn,4))
aknn.model<-paste("KNN",sep="")
aknn<-cbind(aknn.model,aknn)
colnames(aknn)<-colnames(ae_predict)
ae_predict<-rbind(ae_predict,aknn)
pred_knn.model<-paste("KNN",sep="")
pred_knn_class<-cbind(pred_knn.model,t(pred_knn))
colnames(pred_knn_class)<-colnames(class_result)
class_result<-rbind(class_result,pred_knn_class)
#Multiple Paerceptron Network by Stachastic Gradient Descent
nncontrol=trainControl(method = "cv", number = 10) 
  nn_Grid<-expand.grid(
    .size=c(50,5),
    .decay=0.00147)
  nn_fit<-caret::train(label~.,
                       RFtrain,
                       method="nnet",
                       metric="Accuracy",
                       tuneGrid=nn_Grid,
                       MaxNWts=10000,
                       maxit=100,
                       trace=FALSE)
  pred_nn<-predict(nn_fit,test_MM)
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
##Stacked model
#combine all the predictions of above classifiers
label<-as.factor(test_MM$label)
combo<-data.frame(pred,pred_rf,pred_knn,label)
fit_stacked<-caret::train(as.factor(label)~.,
                          combo,
                          method="ranger",
                          tuneGrid=expand.grid(
                       .mtry=1),
                          metric="Accuracy",
                          trControl=fitcontrol)
pred_stacked<-predict(fit_stacked,combo)
#prob_stacked<-predict(fit_stacked,combo,type='prob')
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
pred_stacked<-cbind(pred_stacked.model,t(pred_stacked))
colnames(pred_stacked)<-colnames(class_result)
class_result<-rbind(class_result,pred_stacked)


write.csv(class_result,file="Joaquin_pred_HR_FLAG.csv")
write.csv(ae_predict,file="Joaquin_model_Accuracy.csv")



