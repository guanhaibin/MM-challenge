source("https://bioconductor.org/biocLite.R")
biocLite("Linnorm")
library(Linnorm)
biocLite("limma")
library(limma)
biocLite("biomaRt")
library(biomaRt)
setwd("U:/hguan003/MM")
#Tran_TPM<-read.table("MMRF_CoMMpass_IA9_E74GTF_Salmon_Transcript_TPM.txt",header=T,row.names=1)

Gene_TPM<-read.table("MMRF_CoMMpass_IA9_E74GTF_Salmon_Gene_TPM.txt",header=T,row.names=1)
clinical<-read.csv("globalClinTraining.csv",header=T,sep=",",row.names=2,stringsAsFactors = FALSE)
clinical<-subset(clinical,Study=="MMRF")
clinical<-subset(clinical[,c("Study","HR_FLAG")])
clinical[clinical[,"HR_FLAG"]==TRUE,2]=1
clinical[clinical[,"HR_FLAG"]==FALSE,2]=-1
clin<-clinical[clinical$HR_FLAG!="CENSORED",]

#Convert Ensemble Gene ID to Entez Gene ID 
Transformed<-Linnorm(Gene_TPM)
EnsembleID<-rownames(Transformed)

##biomaRt
ensembl=useMart("ensembl")
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl) 
map=getBM(attributes=c('ensembl_gene_id',"entrezgene"),values="*", mart= ensembl)
EnsembleID<-as.data.frame(EnsembleID)
EnsembleID$entrez<-map$entrezgene[match(EnsembleID$EnsembleID,map$ensembl_gene_id)]
EnsembleID<-EnsembleID[!is.na(EnsembleID$entrez),]
T<-subset(Transformed,rownames(Transformed)%in%EnsembleID[,1])
rownames(T)<-EnsembleID[match(rownames(T),EnsembleID[,1]),2]
colnames(T)<-gsub("_1_BM","",colnames(T))
T<-t(T)
T<-subset(T,rownames(T)%in%rownames(clin))
MMRF<-merge(T,clin,by="row.names",all.x=TRUE)
row.names(MMRF)<-MMRF[,1]
MMRF<-MMRF[-c(1)]
MMRF$Study<-NULL
colnames(MMRF)[ncol(MMRF)]<-c("label")
MMRF$label<-as.factor(MMRF$label)
MMRF<-MMRF[,as.logical(colSums(MMRF != 0))]
write.csv(MMRF,"MMRF_processed.csv")

EnsembleID_to_EntrezID<-EnsembleID
save(EnsembleID_to_EntrezID,file="Ensemble_Map_EntrezID.Rd")
load("Ensemble_Map_EntrezID.Rd")

U_genes_list<-read.csv("U_genes.csv",row.names = 1)
U_genes_list<-as.character(U_genes_list[,1])
E_genes_list<-read.csv("E_genes.csv",row.names = 1)
E_genes_list<-as.character(E_genes_list[,1])
H_genes_list<-read.csv("H_genes.csv",row.names = 1)
H_genes_list<-as.character(H_genes_list[,1])
M_genes_list<-read.csv("M_genes.csv",row.names = 1)
M_genes_list<-as.character(M_genes_list[,1])

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

H_svmfit<-readRDS("./H_svmfit.rds")
H_rffit<-readRDS("./H_rffit.rds")
H_nn_fit<-readRDS("./H_nn_fit.rds")
H_lvq_fit<-readRDS("./H_lvq_fit.rds")
H_gbm_fit<-readRDS("./H_gbm_fit.rds")
H_glmfit<-readRDS("./H_glmfit.rds")
H_fit_stacked<-readRDS("./H_fit_stacked.rds")

M_svmfit<-readRDS("./M_svmfit.rds")
M_rffit<-readRDS("./M_rffit.rds")
M_nn_fit<-readRDS("./M_nn_fit.rds")
M_lvq_fit<-readRDS("./M_lvq_fit.rds")
M_gbm_fit<-readRDS("./M_gbm_fit.rds")
M_glmfit<-readRDS("./M_glmfit.rds")
M_fit_stacked<-readRDS("./M_fit_stacked.rds")

ensemble_rf<-readRDS("./ensemble_rf.rds")


save(U_svmfit,U_rffit,U_nn_fit,U_lvq_fit,U_gbm_fit,U_glmfit,U_fit_stacked,
     E_svmfit,E_rffit,E_nn_fit,E_lvq_fit,E_gbm_fit,E_glmfit,E_fit_stacked,
     H_svmfit,H_rffit,H_nn_fit,H_lvq_fit,H_gbm_fit,H_glmfit,H_fit_stacked,
     M_svmfit,M_rffit,M_nn_fit,M_lvq_fit,M_gbm_fit,M_glmfit,M_fit_stacked,
     U_genes_list,E_genes_list,H_genes_list,M_genes_list,
     ensemble_rf,file="model-state-metadata.Rd")

load("model-state-metadata.Rd")










