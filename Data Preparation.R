#Read data in R and preparation
setwd("U:/hguan003/MM")
clinical<-read.csv("globalClinTraining.csv",header=T,sep=",",row.names = 2,stringsAsFactors = FALSE)
a<-c("GSE24080UAMS","EMTAB4032","HOVON65")
clinical<-subset(clinical,Study%in%a)
clinical<-subset(clinical[,c("Study","D_Age","D_OS","D_OS_FLAG","D_PFS","D_PFS_FLAG","D_ISS","HR_FLAG")])
clinical$Patient<-rownames(clinical)
clinical

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
clinical<-subset(clinical[,c("Patient","Study","D_Age","D_ISS","RS")])
