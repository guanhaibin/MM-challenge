library(glmnet)
library(xgboost)
library(plyr)
library(FCNN4R)
library(MLmetrics)
library(ModelMetrics)
library(stringr)
library(e1071)
library(caret)


setwd("/Users/guanhaibin/Documents/MM")
U_genes_list<-read.csv("U_genes.csv",row.names = 1)
U_genes_list<-as.character(U_genes_list[,1])
E_genes_list<-read.csv("E_genes.csv",row.names = 1)
E_genes_list<-as.character(E_genes_list[,1])
H_genes_list<-read.csv("H_genes.csv",row.names = 1)
H_genes_list<-as.character(H_genes_list[,1])
M_genes_list<-read.csv("M_genes.csv",row.names = 1)
M_genes_list<-as.character(M_genes_list[,1])
load("emc92_uams70.Rd")
genes_list<-unique(c(emc92,uams70))
eu_genes_list<-unique(genes_list)


#read UAMS datasets with selected genes
Gene_expression=read.csv("UAMS_balanced_Ensemble.csv",header=T,sep=",",row.names=1)
Gene_expression<-as.data.frame(t(scale(Gene_expression)))
HR_FLAG<-read.csv("UAMS_target.csv",header=T,sep=",")
HR_FLAG<-t(HR_FLAG)
train_data<-merge(Gene_expression,HR_FLAG,by="row.names",all.x=TRUE)
row.names(train_data)<-train_data[,1]
train_data<-train_data[-c(1)]
colnames(train_data)[ncol(train_data)]<-c("label")
train_data$label<-as.factor(train_data$label)
rownames(train_data)<-gsub("X","U",rownames(train_data))
U_test_MM<-read.csv("U_test_MM.csv",row.names = 1)
U_test_MM<-as.character(U_test_MM[,1])
UAMS<-subset(train_data,rownames(train_data)%in%U_test_MM)

#split UAMS into train ans test set
set.seed(333)
train.rows<-createDataPartition(y=UAMS$label,p=0.7,list=FALSE)
UAMS<-UAMS[train.rows,]
UAMS_test<-UAMS[-train.rows,]


#split UAMS train data into different subset to fit in different models
UAMS_U<-subset(UAMS,select=c(U_genes_list))
UAMS_U$label<-UAMS$label
UAMS_E<-setNames(as.data.frame(matrix(0,nrow=nrow(UAMS),ncol=length(E_genes_list))),E_genes_list)
###match and pick
idx<-colnames(UAMS)[colnames(UAMS)%in%colnames(UAMS_E)]
UAMS_E[,idx]<-UAMS[,idx]
rownames(UAMS_E)<-rownames(UAMS)
UAMS_E$label<-UAMS$label
UAMS_H<-subset(UAMS,select=c(H_genes_list))
UAMS_H$label<-UAMS$label
UAMS_M<-setNames(as.data.frame(matrix(0,nrow=nrow(UAMS),ncol=length(M_genes_list))),M_genes_list)
###match and pick
idx<-colnames(UAMS)[colnames(UAMS)%in%colnames(UAMS_M)]
UAMS_M[,idx]<-UAMS[,idx]
rownames(UAMS_M)<-rownames(UAMS)
UAMS_M$label<-UAMS$label
UAMS_eu<-setNames(as.data.frame(matrix(0,nrow=nrow(UAMS),ncol=length(eu_genes_list))),eu_genes_list)
###match and pick
idx<-colnames(UAMS)[colnames(UAMS)%in%colnames(UAMS_eu)]
UAMS_eu[,idx]<-UAMS[,idx]
rownames(UAMS_eu)<-rownames(UAMS)
UAMS_eu$label<-UAMS$label

#split UAMS test data 
UAMS_U_test<-subset(UAMS_test,select=c(U_genes_list))
UAMS_U_test$label<-UAMS_test$label
UAMS_E_test<-setNames(as.data.frame(matrix(0,nrow=nrow(UAMS_test),ncol=length(E_genes_list))),E_genes_list)
idx<-colnames(UAMS_test)[colnames(UAMS_test)%in%colnames(UAMS_E_test)]
UAMS_E_test[,idx]<-UAMS_test[,idx]
rownames(UAMS_E_test)<-rownames(UAMS_test)
UAMS_E_test$label<-UAMS_test$label
UAMS_H_test<-subset(UAMS_test,select=c(H_genes_list))
UAMS_H_test$label<-UAMS_test$label
UAMS_M_test<-setNames(as.data.frame(matrix(0,nrow=nrow(UAMS_test),ncol=length(M_genes_list))),M_genes_list)
idx<-colnames(UAMS_test)[colnames(UAMS_test)%in%colnames(UAMS_M_test)]
UAMS_M_test[,idx]<-UAMS_test[,idx]
rownames(UAMS_M_test)<-rownames(UAMS_test)
UAMS_M_test$label<-UAMS_test$label
UAMS_eu_test<-setNames(as.data.frame(matrix(0,nrow=nrow(UAMS_test),ncol=length(eu_genes_list))),eu_genes_list)
idx<-colnames(UAMS_test)[colnames(UAMS_test)%in%colnames(UAMS_eu_test)]
UAMS_eu_test[,idx]<-UAMS_test[,idx]
rownames(UAMS_eu_test)<-rownames(UAMS_test)
UAMS_eu_test$label<-UAMS_test$label

#EMTAB
Gene_expression=read.csv("EMTAB_balanced_Ensemble.csv",header=T,sep=",",row.names=1)
Gene_expression<-as.data.frame(t(scale(Gene_expression)))
HR_FLAG<-read.csv("EMTAB_target.csv",header=T,sep=",")
HR_FLAG<-t(HR_FLAG)
train_data<-merge(Gene_expression,HR_FLAG,by="row.names",all.x=TRUE)
row.names(train_data)<-train_data[,1]
train_data<-train_data[-c(1)]
colnames(train_data)[ncol(train_data)]<-c("label")
train_data$label<-as.factor(train_data$label)
rownames(train_data)<-gsub("X","E",rownames(train_data))
E_test_MM<-read.csv("E_test_MM.csv",row.names = 1)
E_test_MM<-as.character(E_test_MM[,1])
EMTAB<-subset(train_data,rownames(train_data)%in%E_test_MM)

#split EMTAB into train ans test set
set.seed(133)
train.rows<-createDataPartition(y=EMTAB$label,p=0.7,list=FALSE)
EMTAB<-EMTAB[train.rows,]
EMTAB_test<-EMTAB[-train.rows,]


#split EMTAB train data into different subset to fit in different models
EMTAB_E<-subset(EMTAB,select=c(E_genes_list))
EMTAB_E$label<-EMTAB$label
EMTAB_U<-setNames(as.data.frame(matrix(0,nrow=nrow(EMTAB),ncol=length(U_genes_list))),U_genes_list)
idx<-colnames(EMTAB)[colnames(EMTAB)%in%colnames(EMTAB_U)]
EMTAB_U[,idx]<-EMTAB[,idx]
rownames(EMTAB_U)<-rownames(EMTAB)
EMTAB_U$label<-EMTAB$label
EMTAB_H<-setNames(as.data.frame(matrix(0,nrow=nrow(EMTAB),ncol=length(H_genes_list))),H_genes_list)
idx<-colnames(EMTAB)[colnames(EMTAB)%in%colnames(EMTAB_H)]
EMTAB_H[,idx]<-EMTAB[,idx]
rownames(EMTAB_H)<-rownames(EMTAB)
EMTAB_H$label<-EMTAB$label
EMTAB_M<-setNames(as.data.frame(matrix(0,nrow=nrow(EMTAB),ncol=length(M_genes_list))),M_genes_list)
idx<-colnames(EMTAB)[colnames(EMTAB)%in%colnames(EMTAB_M)]
EMTAB_M[,idx]<-EMTAB[,idx]
rownames(EMTAB_M)<-rownames(EMTAB)
EMTAB_M$label<-EMTAB$label
EMTAB_eu<-setNames(as.data.frame(matrix(0,nrow=nrow(EMTAB),ncol=length(eu_genes_list))),eu_genes_list)
idx<-colnames(EMTAB)[colnames(EMTAB)%in%colnames(EMTAB_eu)]
EMTAB_eu[,idx]<-EMTAB[,idx]
rownames(EMTAB_eu)<-rownames(EMTAB)
EMTAB_eu$label<-EMTAB$label

#split EMTAB test data 
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
EMTAB_M_test<-setNames(as.data.frame(matrix(0,nrow=nrow(EMTAB_test),ncol=length(M_genes_list))),M_genes_list)
idx<-colnames(EMTAB_test)[colnames(EMTAB_test)%in%colnames(EMTAB_M_test)]
EMTAB_M_test[,idx]<-EMTAB_test[,idx]
rownames(EMTAB_M_test)<-rownames(EMTAB_test)
EMTAB_M_test$label<-EMTAB_test$label
EMTAB_eu_test<-setNames(as.data.frame(matrix(0,nrow=nrow(EMTAB_test),ncol=length(eu_genes_list))),eu_genes_list)
idx<-colnames(EMTAB_test)[colnames(EMTAB_test)%in%colnames(EMTAB_eu_test)]
EMTAB_eu_test[,idx]<-EMTAB_test[,idx]
rownames(EMTAB_eu_test)<-rownames(EMTAB_test)
EMTAB_eu_test$label<-EMTAB_test$label

#read HOVON datasets with selected genes
Gene_expression=read.csv("HOVON_balanced_Ensemble.csv",header=T,sep=",",row.names=1)
Gene_expression<-as.data.frame(t(scale(Gene_expression)))
HR_FLAG<-read.csv("HOVON_target.csv",header=T,sep=",")
HR_FLAG<-t(HR_FLAG)
train_data<-merge(Gene_expression,HR_FLAG,by="row.names",all.x=TRUE)
row.names(train_data)<-train_data[,1]
train_data<-train_data[-c(1)]
colnames(train_data)[ncol(train_data)]<-c("label")
train_data$label<-as.factor(train_data$label)
rownames(train_data)<-gsub("X","H",rownames(train_data))
H_test_MM<-read.csv("H_test_MM.csv",row.names = 1)
H_test_MM<-as.character(H_test_MM[,1])
HOVON<-subset(train_data,rownames(train_data)%in%H_test_MM)

#split HOVON into train ans test set
set.seed(233)
train.rows<-createDataPartition(y=HOVON$label,p=0.7,list=FALSE)
HOVON<-HOVON[train.rows,]
HOVON_test<-HOVON[-train.rows,]


#split HOVON train data into different subset to fit in different models
HOVON_U<-subset(HOVON,select=c(U_genes_list))
HOVON_U$label<-HOVON$label
HOVON_E<-setNames(as.data.frame(matrix(0,nrow=nrow(HOVON),ncol=length(E_genes_list))),E_genes_list)
idx<-colnames(HOVON)[colnames(HOVON)%in%colnames(HOVON_E)]
HOVON_E[,idx]<-HOVON[,idx]
rownames(HOVON_E)<-rownames(HOVON)
HOVON_E$label<-HOVON$label
HOVON_H<-subset(HOVON,select=c(H_genes_list))
HOVON_H$label<-HOVON$label
HOVON_M<-setNames(as.data.frame(matrix(0,nrow=nrow(HOVON),ncol=length(M_genes_list))),M_genes_list)
idx<-colnames(HOVON)[colnames(HOVON)%in%colnames(HOVON_M)]
HOVON_M[,idx]<-HOVON[,idx]
rownames(HOVON_M)<-rownames(HOVON)
HOVON_M$label<-HOVON$label
HOVON_eu<-setNames(as.data.frame(matrix(0,nrow=nrow(HOVON),ncol=length(eu_genes_list))),eu_genes_list)
idx<-colnames(HOVON)[colnames(HOVON)%in%colnames(HOVON_eu)]
HOVON_eu[,idx]<-HOVON[,idx]
rownames(HOVON_eu)<-rownames(HOVON)
HOVON_eu$label<-HOVON$label

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
HOVON_M_test<-setNames(as.data.frame(matrix(0,nrow=nrow(HOVON_test),ncol=length(M_genes_list))),M_genes_list)
idx<-colnames(HOVON_test)[colnames(HOVON_test)%in%colnames(HOVON_M_test)]
HOVON_M_test[,idx]<-HOVON_test[,idx]
rownames(HOVON_M_test)<-rownames(HOVON_test)
HOVON_M_test$label<-HOVON_test$label
HOVON_eu_test<-setNames(as.data.frame(matrix(0,nrow=nrow(HOVON_test),ncol=length(eu_genes_list))),eu_genes_list)
idx<-colnames(HOVON_test)[colnames(HOVON_test)%in%colnames(HOVON_eu_test)]
HOVON_eu_test[,idx]<-HOVON_test[,idx]
rownames(HOVON_eu_test)<-rownames(HOVON_test)
HOVON_eu_test$label<-HOVON_test$label


#read MMRF datasets with selected genes
Gene_expression=read.csv("MMRF_CoMMpass_IA9_E74GTF_Salmon_Gene_TPM.csv",header=TRUE,sep=",",row.names=1)
train_data<-as.data.frame(t(scale(Gene_expression)))
clinical<-read.csv("globalClinTraining.csv",header=TRUE,sep=",",row.names = 2,stringsAsFactors = FALSE)
clinical<-subset(clinical,Study=="MMRF")
clinical<-subset(clinical[,c("Study","HR_FLAG","RNASeq_geneLevelExpFileSamplId")])
clinical[clinical[,"HR_FLAG"]==TRUE,2]=1
clinical[clinical[,"HR_FLAG"]==FALSE,2]=-1
k<-clinical
k<-k[!duplicated(k$RNASeq_geneLevelExpFileSamplId),]
k<-k[!is.na(k$RNASeq_geneLevelExpFileSamplId),]
rownames(k)<-k$RNASeq_geneLevelExpFileSamplId
k$RNASeq_geneLevelExpFileSamplId<-NULL
train_data<-merge(train_data,k,by="row.names",all.x=TRUE)
train_data<-train_data[train_data$HR_FLAG!="CENSORED",]
train_data<-subset(train_data,!is.na(train_data[,1]))
row.names(train_data)<-train_data[,1]
train_data<-train_data[-c(1)]
colnames(train_data)[ncol(train_data)]<-c("label")
train_data$label<-as.factor(train_data$label)
MMRF<-train_data

#split MMRF into train ans test set
set.seed(333)
train.rows<-createDataPartition(y=MMRF$label,p=0.7,list=FALSE)
MMRF<-MMRF[train.rows,]
MMRF_test<-MMRF[-train.rows,]


#split MMRF train data into different subset to fit in different models
MMRF_U<-setNames(as.data.frame(matrix(0,nrow=nrow(MMRF),ncol=length(U_genes_list))),U_genes_list)
###match and pick
idx<-colnames(MMRF)[colnames(MMRF)%in%colnames(MMRF_U)]
MMRF_U[,idx]<-MMRF[,idx]
rownames(MMRF_U)<-rownames(MMRF)
MMRF_U$label<-MMRF$label
MMRF_E<-setNames(as.data.frame(matrix(0,nrow=nrow(MMRF),ncol=length(E_genes_list))),E_genes_list)
###match and pick
idx<-colnames(MMRF)[colnames(MMRF)%in%colnames(MMRF_E)]
MMRF_E[,idx]<-MMRF[,idx]
rownames(MMRF_E)<-rownames(MMRF)
MMRF_E$label<-MMRF$label
MMRF_H<-setNames(as.data.frame(matrix(0,nrow=nrow(MMRF),ncol=length(H_genes_list))),H_genes_list)
###match and pick
idx<-colnames(MMRF)[colnames(MMRF)%in%colnames(MMRF_H)]
MMRF_H[,idx]<-MMRF[,idx]
rownames(MMRF_H)<-rownames(MMRF)
MMRF_H$label<-MMRF$label
MMRF_M<-setNames(as.data.frame(matrix(0,nrow=nrow(MMRF),ncol=length(M_genes_list))),M_genes_list)
###match and pick
idx<-colnames(MMRF)[colnames(MMRF)%in%colnames(MMRF_M)]
MMRF_M[,idx]<-MMRF[,idx]
rownames(MMRF_M)<-rownames(MMRF)
MMRF_M$label<-MMRF$label
MMRF_eu<-setNames(as.data.frame(matrix(0,nrow=nrow(MMRF),ncol=length(eu_genes_list))),eu_genes_list)
###match and pick
idx<-colnames(MMRF)[colnames(MMRF)%in%colnames(MMRF_eu)]
MMRF_eu[,idx]<-MMRF[,idx]
rownames(MMRF_eu)<-rownames(MMRF)
MMRF_eu$label<-MMRF$label

#split MMRF test data 
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
MMRF_M_test<-setNames(as.data.frame(matrix(0,nrow=nrow(MMRF_test),ncol=length(M_genes_list))),M_genes_list)
idx<-colnames(MMRF_test)[colnames(MMRF_test)%in%colnames(MMRF_M_test)]
MMRF_M_test[,idx]<-MMRF_test[,idx]
rownames(MMRF_M_test)<-rownames(MMRF_test)
MMRF_M_test$label<-MMRF_test$label
MMRF_eu_test<-setNames(as.data.frame(matrix(0,nrow=nrow(MMRF_test),ncol=length(eu_genes_list))),eu_genes_list)
idx<-colnames(MMRF_test)[colnames(MMRF_test)%in%colnames(MMRF_eu_test)]
MMRF_eu_test[,idx]<-MMRF_test[,idx]
rownames(MMRF_eu_test)<-rownames(MMRF_test)
MMRF_eu_test$label<-MMRF_test$label




##combine all train data corresponding to the different model
train_U<-rbind(UAMS_U,EMTAB_U,HOVON_U,MMRF_U)
train_E<-rbind(UAMS_E,EMTAB_E,HOVON_E,MMRF_E)
train_H<-rbind(UAMS_H,EMTAB_H,HOVON_H,MMRF_H)
train_M<-rbind(UAMS_M,EMTAB_M,HOVON_M,MMRF_M)
train_eu<-rbind(UAMS_eu,EMTAB_eu,HOVON_eu,MMRF_eu)

train_U_label<-train_U$label
train_U$label<-NULL
train_E_label<-train_E$label
train_E$label<-NULL
train_H_label<-train_H$label
train_H$label<-NULL
train_M_label<-train_M$label
train_M$label<-NULL
train_eu_label<-train_eu$label
train_eu$label<-NULL

##combine all test data correspomding to the different model
test_U<-rbind(UAMS_U_test,EMTAB_U_test,HOVON_U_test,MMRF_U_test)
test_E<-rbind(UAMS_E_test,EMTAB_E_test,HOVON_E_test,MMRF_E_test)
test_H<-rbind(UAMS_H_test,EMTAB_H_test,HOVON_H_test,MMRF_H_test)
test_M<-rbind(UAMS_M_test,EMTAB_M_test,HOVON_M_test,MMRF_M_test)
test_eu<-rbind(UAMS_eu_test,EMTAB_eu_test,HOVON_eu_test,MMRF_eu_test)

test_U_label<-test_U$label
test_U$label<-NULL
test_E_label<-test_E$label
test_E$label<-NULL
test_H_label<-test_H$label
test_H$label<-NULL
test_M_label<-test_M$label
test_M$label<-NULL
test_eu_label<-test_eu$label
test_eu$label<-NULL


##impute
myimpute <- function(data,clus = 20) # clus = Number or row in splited matrix
{
  library(impute)
  data <- data.frame(data)
  row_data <-nrow(data)
  ind<-as.factor(c(gl(round(row_data/clus)-1,clus),rep(round(row_data/clus)-1+1, nrow(data)-length(gl(round(row_data/clus)-1,clus)))))
  newMat <- split(data, ind)
  res <- lapply(newMat,function(x)impute.knn(as.matrix(x)))
  res <- lapply(res,"[[","data")
  res <- do.call(rbind, res)
  res
}

a<-myimpute(t(train_U),clus=20)
train_U<-as.data.frame(t(a))
train_U$label<-train_U_label
a<-myimpute(t(test_U),clus=20)
test_U<-as.data.frame(t(a))
test_U$label<-test_U_label
a<-myimpute(t(train_E),clus=20)
train_E<-as.data.frame(t(a))
train_E$label<-train_E_label
a<-myimpute(t(test_E),clus=20)
test_E<-as.data.frame(t(a))
test_E$label<-test_E_label
a<-myimpute(t(train_H),clus=20)
train_H<-as.data.frame(t(a))
train_H$label<-train_H_label
a<-myimpute(t(test_H),clus=20)
test_H<-as.data.frame(t(a))
test_H$label<-test_H_label

a<-myimpute(t(train_M),clus=20)
train_M<-as.data.frame(t(a))
train_M$label<-train_M_label
a<-myimpute(t(test_M),clus=20)
test_M<-as.data.frame(t(a))
test_M$label<-test_M_label
a<-myimpute(t(train_eu),clus=20)
train_eu<-as.data.frame(t(a))
train_eu$label<-train_eu_label
a<-myimpute(t(test_eu),clus=20)
test_eu<-as.data.frame(t(a))
test_eu$label<-test_eu_label



th<-0.5




##load the U_mdoel 
U_svmfit<-readRDS("./U_svmfit.rds")
train_U_svm<-predict(U_svmfit,train_U[,1:(ncol(train_U)-1)])
test_U_svm<-predict(U_svmfit,test_U[,1:(ncol(test_U)-1)])
U_rffit<-readRDS("./U_rffit.rds")
train_U_rf_prob<-predict(U_rffit,train_U[,1:(ncol(train_U)-1)],type="prob")
train_U_rf<-factor(ifelse(train_U_rf_prob$X1>th,"1","-1"))
test_U_rf_prob<-predict(U_rffit,test_U[,1:(ncol(test_U)-1)],type="prob")
test_U_rf<-factor(ifelse(test_U_rf_prob$X1>th,"1","-1"))
U_nn_fit<-readRDS("./U_nn_fit.rds")
train_U_nn_prob<-predict(U_nn_fit,train_U[,1:(ncol(train_U)-1)],type="prob")
train_U_nn<-predict(U_nn_fit,train_U[,1:(ncol(train_U)-1)])
test_U_nn_prob<-predict(U_nn_fit,test_U[,1:(ncol(test_U)-1)],type="prob")
test_U_nn<-predict(U_nn_fit,test_U[,1:(ncol(test_U)-1)])
U_lvq_fit<-readRDS("./U_lvq_fit.rds")
train_U_lvq<-predict(U_lvq_fit,train_U[,1:(ncol(train_U)-1)])
test_U_lvq<-predict(U_lvq_fit,test_U[,1:(ncol(test_U)-1)])
U_gbm_fit<-readRDS("./U_gbm_fit.rds")
train_U_gbm_prob<-predict(U_gbm_fit,train_U[,1:(ncol(train_U)-1)],type="prob")
train_U_gbm<-factor(ifelse(train_U_gbm_prob$X1>th,"1","-1"))
test_U_gbm_prob<-predict(U_gbm_fit,test_U[,1:(ncol(test_U)-1)],type="prob")
test_U_gbm<-factor(ifelse(test_U_gbm_prob$X1>th,"1","-1"))
U_glmfit<-readRDS("./U_glmfit.rds")
train_U_glm_prob<-predict(U_glmfit,train_U[,1:(ncol(train_U)-1)],type="prob")
train_U_glm<-factor(ifelse(train_U_glm_prob$X1>th,"1","-1"))
test_U_glm_prob<-predict(U_glmfit,test_U[,1:(ncol(test_U)-1)],type="prob")
test_U_glm<-factor(ifelse(test_U_glm_prob$X1>th,"1","-1"))
combo_U_train<-data.frame(train_U_svm,train_U_rf,train_U_nn,train_U_lvq,train_U_gbm,train_U_glm,label<-as.factor(train_U$label))
colnames(combo_U_train)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train","train_label" )
U_fit_stacked<-readRDS("./U_fit_stacked.rds")
train_U_prob<-predict(U_fit_stacked,combo_U_train,type="prob")
th<-0.5
train_U_stacked<-factor(ifelse(train_U_prob$X1>th,"1","-1"))
combo_U_test<-data.frame(test_U_svm,test_U_rf,test_U_nn,test_U_lvq,test_U_gbm,test_U_glm,label<-as.factor(test_U$label))
colnames(combo_U_test)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train", "train_label" )
test_U_prob<-predict(U_fit_stacked,combo_U_test,type="prob")
th<-0.5
test_U_stacked<-factor(ifelse(test_U_prob$X1>th,"1","-1"))


##load the E_mdoel
E_svmfit<-readRDS("./E_svmfit.rds")
train_E_svm<-predict(E_svmfit,train_E[,1:(ncol(train_E)-1)])
test_E_svm<-predict(E_svmfit,test_E[,1:(ncol(test_E)-1)])
E_rffit<-readRDS("./E_rffit.rds")
train_E_rf_prob<-predict(E_rffit,train_E[,1:(ncol(train_E)-1)],type="prob")
train_E_rf<-factor(ifelse(train_E_rf_prob$X1>th,"1","-1"))
test_E_rf_prob<-predict(E_rffit,test_E[,1:(ncol(test_E)-1)],type="prob")
test_E_rf<-factor(ifelse(test_E_rf_prob$X1>th,"1","-1"))
E_nn_fit<-readRDS("./E_nn_fit.rds")
train_E_nn_prob<-predict(E_nn_fit,train_E[,1:(ncol(train_E)-1)],type="prob")
train_E_nn<-predict(E_nn_fit,train_E[,1:(ncol(train_E)-1)])
test_E_nn_prob<-predict(E_nn_fit,test_E[,1:(ncol(test_E)-1)],type="prob")
test_E_nn<-predict(E_nn_fit,test_E[,1:(ncol(test_E)-1)])
E_lvq_fit<-readRDS("./E_lvq_fit.rds")
train_E_lvq<-predict(E_lvq_fit,train_E[,1:(ncol(train_E)-1)])
test_E_lvq<-predict(E_lvq_fit,test_E[,1:(ncol(test_E)-1)])
E_gbm_fit<-readRDS("./E_gbm_fit.rds")
train_E_gbm_prob<-predict(E_gbm_fit,train_E[,1:(ncol(train_E)-1)],type="prob")
train_E_gbm<-factor(ifelse(train_E_gbm_prob$X1>th,"1","-1"))
test_E_gbm_prob<-predict(E_gbm_fit,test_E[,1:(ncol(test_E)-1)],type="prob")
test_E_gbm<-factor(ifelse(test_E_gbm_prob$X1>th,"1","-1"))
E_glmfit<-readRDS("./E_glmfit.rds")
train_E_glm_prob<-predict(E_glmfit,train_E[,1:(ncol(train_E)-1)],type="prob")
train_E_glm<-factor(ifelse(train_E_glm_prob$X1>th,"1","-1"))
test_E_glm_prob<-predict(E_glmfit,test_E[,1:(ncol(test_E)-1)],type="prob")
test_E_glm<-factor(ifelse(test_E_glm_prob$X1>th,"1","-1"))
combo_E_train<-data.frame(train_E_svm,train_E_rf,train_E_nn,train_E_lvq,train_E_gbm,train_E_glm,label<-as.factor(train_E$label))
colnames(combo_E_train)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train","train_label" )
E_fit_stacked<-readRDS("./E_fit_stacked.rds")
train_E_prob<-predict(E_fit_stacked,combo_E_train,type="prob")
th<-0.5
train_E_stacked<-factor(ifelse(train_E_prob$X1>th,"1","-1"))
combo_E_test<-data.frame(test_E_svm,test_E_rf,test_E_nn,test_E_lvq,test_E_gbm,test_E_glm,label<-as.factor(test_E$label))
colnames(combo_E_test)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train","train_label" )
test_E_prob<-predict(E_fit_stacked,combo_E_test,type="prob")
th<-0.5
test_E_stacked<-factor(ifelse(test_E_prob$X1>th,"1","-1"))

##load the H_mdoel
H_svmfit<-readRDS("./H_svmfit.rds")
train_H_svm<-predict(H_svmfit,train_H[,1:(ncol(train_H)-1)])
test_H_svm<-predict(H_svmfit,test_H[,1:(ncol(test_H)-1)])
H_rffit<-readRDS("./H_rffit.rds")
train_H_rf_prob<-predict(H_rffit,train_H[,1:(ncol(train_H)-1)],type="prob")
train_H_rf<-factor(ifelse(train_H_rf_prob$X1>th,"1","-1"))
test_H_rf_prob<-predict(H_rffit,test_H[,1:(ncol(test_H)-1)],type="prob")
test_H_rf<-factor(ifelse(test_H_rf_prob$X1>th,"1","-1"))
H_nn_fit<-readRDS("./H_nn_fit.rds")
train_H_nn_prob<-predict(H_nn_fit,train_H[,1:(ncol(train_H)-1)],type="prob")
train_H_nn<-predict(H_nn_fit,train_H[,1:(ncol(train_H)-1)])
test_H_nn_prob<-predict(H_nn_fit,test_H[,1:(ncol(test_H)-1)],type="prob")
test_H_nn<-predict(H_nn_fit,test_H[,1:(ncol(test_H)-1)])
H_lvq_fit<-readRDS("./H_lvq_fit.rds")
train_H_lvq<-predict(H_lvq_fit,train_H[,1:(ncol(train_H)-1)])
test_H_lvq<-predict(H_lvq_fit,test_H[,1:(ncol(test_H)-1)])
H_gbm_fit<-readRDS("./H_gbm_fit.rds")
train_H_gbm_prob<-predict(H_gbm_fit,train_H[,1:(ncol(train_H)-1)],type="prob")
train_H_gbm<-factor(ifelse(train_H_gbm_prob$X1>th,"1","-1"))
test_H_gbm_prob<-predict(H_gbm_fit,test_H[,1:(ncol(test_H)-1)],type="prob")
test_H_gbm<-factor(ifelse(test_H_gbm_prob$X1>th,"1","-1"))
H_glmfit<-readRDS("./H_glmfit.rds")
train_H_glm_prob<-predict(H_glmfit,train_H[,1:(ncol(train_H)-1)],type="prob")
train_H_glm<-factor(ifelse(train_H_glm_prob$X1>th,"1","-1"))
test_H_glm_prob<-predict(H_glmfit,test_H[,1:(ncol(test_H)-1)],type="prob")
test_H_glm<-factor(ifelse(test_H_glm_prob$X1>th,"1","-1"))
combo_H_train<-data.frame(train_H_svm,train_H_rf,train_H_nn,train_H_lvq,train_H_gbm,train_H_glm,label<-as.factor(train_H$label))
colnames(combo_H_train)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train","train_label" )
H_fit_stacked<-readRDS("./H_fit_stacked.rds")
train_H_prob<-predict(H_fit_stacked,combo_H_train,type="prob")
th<-0.5
train_H_stacked<-factor(ifelse(train_H_prob$X1>th,"1","-1"))
combo_H_test<-data.frame(test_H_svm,test_H_rf,test_H_nn,test_H_lvq,test_H_gbm,test_H_glm,label<-as.factor(test_H$label))
colnames(combo_H_test)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train","train_label" )
test_H_prob<-predict(H_fit_stacked,combo_H_test,type="prob")
th<-0.5
test_H_stacked<-factor(ifelse(test_H_prob$X1>th,"1","-1"))


##load the M_mdoel
M_svmfit<-readRDS("./M_svmfit.rds")
train_M_svm<-predict(M_svmfit,train_M[,1:(ncol(train_M)-1)])
test_M_svm<-predict(M_svmfit,test_M[,1:(ncol(test_M)-1)])
M_rffit<-readRDS("./M_rffit.rds")
train_M_rf_prob<-predict(M_rffit,train_M[,1:(ncol(train_M)-1)],type="prob")
train_M_rf<-factor(ifelse(train_M_rf_prob$X1>th,"1","-1"))
test_M_rf_prob<-predict(M_rffit,test_M[,1:(ncol(test_M)-1)],type="prob")
test_M_rf<-factor(ifelse(test_M_rf_prob$X1>th,"1","-1"))
M_nn_fit<-readRDS("./M_nn_fit.rds")
train_M_nn_prob<-predict(M_nn_fit,train_M[,1:(ncol(train_M)-1)],type="prob")
train_M_nn<-predict(M_nn_fit,train_M[,1:(ncol(train_M)-1)])
test_M_nn_prob<-predict(M_nn_fit,test_M[,1:(ncol(test_M)-1)],type="prob")
test_M_nn<-predict(M_nn_fit,test_M[,1:(ncol(test_M)-1)])
M_lvq_fit<-readRDS("./M_lvq_fit.rds")
train_M_lvq<-predict(M_lvq_fit,train_M[,1:(ncol(train_M)-1)])
test_M_lvq<-predict(M_lvq_fit,test_M[,1:(ncol(test_M)-1)])
M_gbm_fit<-readRDS("./M_gbm_fit.rds")
train_M_gbm_prob<-predict(M_gbm_fit,train_M[,1:(ncol(train_M)-1)],type="prob")
train_M_gbm<-factor(ifelse(train_M_gbm_prob$X1>th,"1","-1"))
test_M_gbm_prob<-predict(M_gbm_fit,test_M[,1:(ncol(test_M)-1)],type="prob")
test_M_gbm<-factor(ifelse(test_M_gbm_prob$X1>th,"1","-1"))
M_glmfit<-readRDS("./M_glmfit.rds")
train_M_glm_prob<-predict(M_glmfit,train_M[,1:(ncol(train_M)-1)],type="prob")
train_M_glm<-factor(ifelse(train_M_glm_prob$X1>th,"1","-1"))
test_M_glm_prob<-predict(M_glmfit,test_M[,1:(ncol(test_M)-1)],type="prob")
test_M_glm<-factor(ifelse(test_M_glm_prob$X1>th,"1","-1"))
combo_M_train<-data.frame(train_M_svm,train_M_rf,train_M_nn,train_M_lvq,train_M_gbm,train_M_glm,label<-as.factor(train_M$label))
colnames(combo_M_train)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train","train_label" )
M_fit_stacked<-readRDS("./M_fit_stacked.rds")
train_M_prob<-predict(M_fit_stacked,combo_M_train,type="prob")
th<-0.5
train_M_stacked<-factor(ifelse(train_M_prob$X1>th,"1","-1"))
combo_M_test<-data.frame(test_M_svm,test_M_rf,test_M_nn,test_M_lvq,test_M_gbm,test_M_glm,label<-as.factor(test_M$label))
colnames(combo_M_test)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train","train_label" )
test_M_prob<-predict(M_fit_stacked,combo_M_test,type="prob")
th<-0.5
test_M_stacked<-factor(ifelse(test_M_prob$X1>th,"1","-1"))





##load the eu_mdoel 
eu_svmfit<-readRDS("./eu_svmfit.rds")
train_eu_svm<-predict(eu_svmfit,train_eu[,1:(ncol(train_eu)-1)])
test_eu_svm<-predict(eu_svmfit,test_eu[,1:(ncol(test_eu)-1)])
eu_rffit<-readRDS("./eu_rffit.rds")
train_eu_rf_prob<-predict(eu_rffit,train_eu[,1:(ncol(train_eu)-1)],type="prob")
train_eu_rf<-factor(ifelse(train_eu_rf_prob$X1>th,"1","-1"))
test_eu_rf_prob<-predict(eu_rffit,test_eu[,1:(ncol(test_eu)-1)],type="prob")
test_eu_rf<-factor(ifelse(test_eu_rf_prob$X1>th,"1","-1"))
eu_nn_fit<-readRDS("./eu_nn_fit.rds")
train_eu_nn_prob<-predict(eu_nn_fit,train_eu[,1:(ncol(train_eu)-1)],type="prob")
train_eu_nn<-predict(eu_nn_fit,train_eu[,1:(ncol(train_eu)-1)])
test_eu_nn_prob<-predict(eu_nn_fit,test_eu[,1:(ncol(test_eu)-1)],type="prob")
test_eu_nn<-predict(eu_nn_fit,test_eu[,1:(ncol(test_eu)-1)])
eu_lvq_fit<-readRDS("./eu_lvq_fit.rds")
train_eu_lvq<-predict(eu_lvq_fit,train_eu[,1:(ncol(train_eu)-1)])
test_eu_lvq<-predict(eu_lvq_fit,test_eu[,1:(ncol(test_eu)-1)])
eu_gbm_fit<-readRDS("./eu_gbm_fit.rds")
train_eu_gbm_prob<-predict(eu_gbm_fit,train_eu[,1:(ncol(train_eu)-1)],type="prob")
train_eu_gbm<-factor(ifelse(train_eu_gbm_prob$X1>th,"1","-1"))
test_eu_gbm_prob<-predict(eu_gbm_fit,test_eu[,1:(ncol(test_eu)-1)],type="prob")
test_eu_gbm<-factor(ifelse(test_eu_gbm_prob$X1>th,"1","-1"))
eu_glmfit<-readRDS("./eu_glmfit.rds")
train_eu_glm_prob<-predict(eu_glmfit,train_eu[,1:(ncol(train_eu)-1)],type="prob")
train_eu_glm<-factor(ifelse(train_eu_glm_prob$X1>th,"1","-1"))
test_eu_glm_prob<-predict(eu_glmfit,test_eu[,1:(ncol(test_eu)-1)],type="prob")
test_eu_glm<-factor(ifelse(test_eu_glm_prob$X1>th,"1","-1"))
combo_eu_train<-data.frame(train_eu_svm,train_eu_rf,train_eu_nn,train_eu_lvq,train_eu_gbm,train_eu_glm,label<-as.factor(train_eu$label))
colnames(combo_eu_train)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train","train_label" )
eu_fit_stacked<-readRDS("./eu_fit_stacked.rds")
train_eu_prob<-predict(eu_fit_stacked,combo_eu_train,type="prob")
th<-0.5
train_eu_stacked<-factor(ifelse(train_eu_prob$X1>th,"1","-1"))
combo_eu_test<-data.frame(test_eu_svm,test_eu_rf,test_eu_nn,test_eu_lvq,test_eu_gbm,test_eu_glm,label<-as.factor(test_eu$label))
colnames(combo_eu_test)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train", "train_label" )
test_eu_prob<-predict(eu_fit_stacked,combo_eu_test,type="prob")
th<-0.5
test_eu_stacked<-factor(ifelse(test_eu_prob$X1>th,"1","-1"))




##global ensemble_train data
U_models_train<-cbind(combo_U_train[,1:6],train_U_stacked)
colnames(U_models_train)<-gsub("_pred_","_U_",colnames(U_models_train))
E_models_train<-cbind(combo_E_train[,1:6],train_E_stacked)
colnames(E_models_train)<-gsub("_pred_","_E_",colnames(E_models_train))
H_models_train<-cbind(combo_H_train[,1:6],train_H_stacked)
colnames(H_models_train)<-gsub("_pred_","_H_",colnames(H_models_train))
M_models_train<-cbind(combo_M_train[,1:6],train_M_stacked)
colnames(M_models_train)<-gsub("_pred_","_M_",colnames(M_models_train))
eu_models_train<-cbind(combo_eu_train[,1:6],train_eu_stacked)
colnames(eu_models_train)<-gsub("_pred_","_eu_",colnames(eu_models_train))

ensemble_train<-cbind(U_models_train,E_models_train,H_models_train,M_models_train,eu_models_train)
ensemble_train$label<-combo_U_train$train_label


##global ensemble_test data
U_models_test<-cbind(combo_U_test[,1:6],test_U_stacked)
colnames(U_models_test)<-gsub("_pred_","_U_",colnames(U_models_test))
E_models_test<-cbind(combo_E_test[,1:6],test_E_stacked)
colnames(E_models_test)<-gsub("_pred_","_E_",colnames(E_models_test))
H_models_test<-cbind(combo_H_test[,1:6],test_H_stacked)
colnames(H_models_test)<-gsub("_pred_","_H_",colnames(H_models_test))
M_models_test<-cbind(combo_M_test[,1:6],test_M_stacked)
colnames(M_models_test)<-gsub("_pred_","_M_",colnames(M_models_test))
eu_models_test<-cbind(combo_eu_test[,1:6],test_eu_stacked)
colnames(eu_models_test)<-gsub("_pred_","_eu_",colnames(eu_models_test))

ensemble_test<-cbind(U_models_test,E_models_test,H_models_test,M_models_test,eu_models_test)
ensemble_test$label<-combo_U_test$train_label
colnames(ensemble_test)<-colnames(ensemble_train)

##train global ensemble model
#Random Forest
set.seed(7)
fitcontrol<-trainControl(method="cv",number=10,classProbs = TRUE) ###
levels(ensemble_train$label)<-make.names(levels(factor(ensemble_train$label)))
for (i in 1:(ncol(ensemble_train)-1)){
levels(ensemble_train[,i])<-c("-1","1")
}
ensemble_rf<-caret::train(as.factor(label)~.,
                       ensemble_train,
                       method="parRF",
                       tunGrid=expand.grid(
                         .mtry=1),
                       metric="Accuracy",
                       trControl=fitcontrol)
saveRDS(ensemble_rf,"./ensemble_rf.rds")
ensemble_rf_prob<-predict(ensemble_rf,ensemble_test[,1:35],type="prob")  ###
th<-0.5
ensemble_rf_pred<-factor(ifelse(ensemble_rf_prob$X1>th,"1","-1"))
ensemble_rf_accuracy<-Accuracy(ensemble_rf_pred,ensemble_test$label)

#Neural Network
nn_Grid<-expand.grid(
  .size=c(10,5,3,2),
  .decay=0.00147)
ensemble_nn<-caret::train(as.factor(label)~.,
                          ensemble_train,
                          method="nnet",
                          tuneGrid=nn_Grid,
                          metric="Accuracy",
                          MaxNWts=10000,
                          maxit=100,
                          trace=FALSE,
                          trControl=fitcontrol)
##save ensemble_rf model 
saveRDS(ensemble_nn,"./ensemble_nn.rds")
ensemble_nn_prob<-predict(ensemble_nn,ensemble_test[,1:35],type="prob")
th<-0.5
ensemble_nn_pred<-factor(ifelse(ensemble_nn_prob$X1>th,"1","-1"))
ensemble_nn_accuracy<-Accuracy(ensemble_nn_pred,ensemble_test$label)


ensemble_accuracy<-as.data.frame(matrix(0,ncol=1,nrow = 37))
for (i in 1:(ncol(ensemble_test)-1)){
  levels(ensemble_test[,i])<-levels(ensemble_test$label)
  a<-Accuracy(ensemble_test[,i],ensemble_test$label)
  ensemble_accuracy[i,1]<-a
  rownames(ensemble_accuracy)[i]<-colnames(ensemble_test)[i]
}
ensemble_accuracy[36,1]<-ensemble_rf_accuracy
rownames(ensemble_accuracy)[36]<-"ensemble_rf"
ensemble_accuracy[37,1]<-ensemble_nn_accuracy
rownames(ensemble_accuracy)[37]<-"ensemble_nn"
rownames(ensemble_accuracy)<-gsub("train","test",rownames(ensemble_accuracy))
write.csv(ensemble_accuracy,"ensemble_accuracy.csv")







