

#read UAMS datasets with selected genes
setwd("U:/hguan003/MM")
Gene_expression=read.csv("GSE24080UAMSentrezIDlevel.csv",header=T,sep=",",row.names=1)
clinical<-read.csv("globalClinTraining.csv",header=T,sep=",",row.names = 2,stringsAsFactors = FALSE)
a<-c("GSE24080UAMS","EMTAB4032","HOVON65")
clinical<-subset(clinical,Study%in%a)
clinical<-subset(clinical[,c("Study","D_Age","D_PFS","D_PFS_FLAG","D_ISS","HR_FLAG")])
clinical$Patient<-rownames(clinical)
clinical<-clinical[c("D_Age","D_ISS","Study","Patient","D_PFS","D_PFS_FLAG","HR_FLAG")]
clinical$D_Age<-clinical$D_Age/100
clinical$D_ISS<-clinical$D_ISS/10
clinical_UAMS<-subset(clinical,Study=="GSE24080UAMS")
clinical_UAMS[clinical_UAMS[,"HR_FLAG"]==TRUE,7]=1
clinical_UAMS[clinical_UAMS[,"HR_FLAG"]==FALSE,7]=0
train_data=as.data.frame(t(Gene_expression))
train_data<-as.data.frame(scale(train_data))
train_data<-merge(train_data,clinical_UAMS[,c("D_Age","D_ISS","Study","Patient","D_PFS","D_PFS_FLAG","HR_FLAG")],by="row.names",all.x=TRUE)
row.names(train_data)<-train_data[,1]
train_data<-train_data[-c(1)]

