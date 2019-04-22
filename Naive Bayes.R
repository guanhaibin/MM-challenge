#Naive Bayes
set.seed(7)
nb_fit<-caret::train(label~.,
                     RFtrain,
                     method="nb",
                     metric="Accuracy",
                     trControl=fitcontrol)
pred_nb<-predict(nb_fit,test_MM[,1:nfea])
nb_pred_train<-predict(nb_fit,RFtrain[,1:nfea])
pred_nb<-as.vector(pred_nb)
test_MM$label<-as.vector(test_MM$label)
nb_accuracy<-Accuracy(pred_nb,test_MM$label)
nb_AUC<-adjustAUC(pred_nb,as.numeric(test_MM$label))
nb_F1<-F1_Score(as.numeric(test_MM$label),as.numeric(pred_nb))
nb_PRAUC<-PRAUC(as.numeric(pred_nb),as.numeric(test_MM$label))
nb_MCC<-adjustmcc(as.numeric(test_MM$label),as.numeric(pred_nb),cutoff=1)
nb_sensitivity<-Sensitivity(test_MM$label,pred_nb)
nb_specific<-Specificity(test_MM$label,pred_nb)
nb_classloss<-ZeroOneLoss(test_MM$label,pred_nb)
#add model gbm prediction score to class_result and ae_predict
a_nb<-as.data.frame(t(c(nb_accuracy,nb_AUC,nb_F1,nb_PRAUC,nb_MCC,nb_sensitivity,nb_specific,nb_classloss)))
a_nb<-as.data.frame(round(a_nb,4))
a_nb.model<-paste("NB",sep="")
a_nb<-cbind(a_nb.model,a_nb)
colnames(a_nb)<-colnames(ae_predict)
ae_predict<-rbind(ae_predict,a_nb)
pred_nb.model<-paste("NB",sep="")
pred_nb_class<-cbind(pred_nb.model,t(pred_nb))
colnames(pred_nb_class)<-colnames(class_result)
class_result<-rbind(class_result,pred_nb_class)


U_nb_fit<-readRDS("./U_nb_fit.rds")
train_U_nb<-predict(U_nb_fit,train_U[,1:(ncol(train_U)-1)])
test_U_nb<-predict(U_nb_fit,test_U[,1:(ncol(test_U)-1)])

E_nb_fit<-readRDS("./E_nb_fit.rds")
train_E_nb<-predict(E_nb_fit,train_E[,1:(ncol(train_E)-1)])
test_E_nb<-predict(E_nb_fit,test_E[,1:(ncol(test_E)-1)])

H_nb_fit<-readRDS("./H_nb_fit.rds")
train_H_nb<-predict(H_nb_fit,train_H[,1:(ncol(train_H)-1)])
test_H_nb<-predict(H_nb_fit,test_H[,1:(ncol(test_H)-1)])


U_nb_fit<-readRDS("./U_nb_fit.rds")
test_U_nb<-predict(U_nb_fit,test_U[,1:(ncol(test_U)-1)])

E_nb_fit<-readRDS("./E_nb_fit.rds")
test_E_nb<-predict(E_nb_fit,test_E[,1:(ncol(test_E)-1)])

H_nb_fit<-readRDS("./H_nb_fit.rds")
test_H_nb<-predict(H_nb_fit,test_H[,1:(ncol(test_H)-1)])

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







#####match specific genes

for (i in 1:ncol(UAMS_E_test)){
  for (j in 1:ncol(UAMS_test)) {
    if (colnames(UAMS_test)[j]==colnames(UAMS_E_test)[i]) {
      UAMS_E_test[,i]<-UAMS_test[,j]
    }
  }
}






##eXtreme Gradient Boosting package
install.packages("xgboost")
library(xgboost)
install.packages("Matrix")
library(Matrix)



colnames(train_MM)[1:nfea]<-paste("g",colnames(train_MM)[1:nfea],sep='')
train_MM$label<-as.numeric(train_MM$label)
for (i in 1:length(train_MM$label)) {if (train_MM$label[i]==-1 ){train_MM$label[i]=0}}
train_MM$label<-as.factor(train_MM$label)
xg_train<-sparse.model.matrix(~label,data=train_MM)
colnames(test_MM)<-paste("g",colnames(test_MM),sep='')
test_MM$label<-1
xg_test<-sparse.model.matrix(label~.-1,data=test_MM)
###XGBoost
xg.model <- xgboost(data= xg_train,     #train sparse matrix 
                    label= train_MM$label,          #output vector to be predicted 
                    eval.metric= 'logloss',        #model minimizes Root Mean Squared Error
                    objective= "binary:logistic",     #regression
                    #tuning parameters
                    max.depth= 8,            #Vary btwn 3-15
                    eta= 0.1,                #Vary btwn 0.1-0.3
                    nthread = 5,             #Increase this to improve speed
                    subsample= 1,            #Vary btwn 0.8-1
                    colsample_bytree= 0.5,   #Vary btwn 0.3-0.8
                    lambda= 0.5,             #Vary between 0-3
                    alpha= 0.5,              #Vary between 0-3
                    min_child_weight= 3,     #Vary btwn 1-10
                    nround= 30 )             

pred_xg<-predict(xg.model,xg_test)






setwd("//bear-ad.cs.fiu.edu/homes/Downloads")
a<-load("model-state-metadata.Rd")
str(model.state)
class(model.state)







