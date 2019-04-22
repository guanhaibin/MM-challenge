#!/usr/bin/env Rscript

##R library dependencies
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("Matrix.utils"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("caret"))
suppressPackageStartupMessages(library("e1071"))
suppressPackageStartupMessages(library("gbm"))
suppressPackageStartupMessages(library("class"))
suppressPackageStartupMessages(library("glmnet"))
suppressPackageStartupMessages(library("ranger"))
suppressPackageStartupMessages(library("randomForest"))
suppressPackageStartupMessages(library("FCNN4R"))
suppressPackageStartupMessages(library("Linnorm"))
suppressPackageStartupMessages(library("xgboost"))
suppressPackageStartupMessages(library("import"))



use.sequestered.data <- FALSE

if(use.sequestered.data) {
  suppressPackageStartupMessages(library("synapseClient"))
  synapseLogin()
}


## This R script assumes that it runs in the root directory of the Docker image and that
## the test data are mounted in ./test-data,
## the output should be written to ./output,
## and that the entire directory structure of the submitted Docker image
## (e.g., an R object encapsulating trained modeler state) is mounted at ./
docker.image.base <- "./"
test.dir <- "./test-data/"
output.dir <- "./output/"


## Read in the trained model state.
model.state.metadata.file <- paste0(docker.image.base, "model-state-metadata.Rd")
load(model.state.metadata.file)


## This is a kludge to run this without Docker.  Download only the files we are
## interested in from the sequester site (at synId = syn9763945)
chal_data_df <- NULL

if(use.sequestered.data) {
  chal_data_table <- synTableQuery('select id,name from syn9763945')
  chal_data_df <- chal_data_table@values
  
  ## Download the annotation file
  chal_data_df_anno <- subset(chal_data_df, grepl(name, pattern="sc2_Validation_ClinAnnotations.csv"))
  sapply(chal_data_df_anno$id, function(x) { synGet(x, downloadLocation=test.dir)})
}


## Read in the annotation file
training.validation.file <- paste0(test.dir, "sc2_Validation_ClinAnnotations.csv")
anno.tbl <- read.table(training.validation.file, sep=",", header=TRUE,row.names = 2,stringsAsFactors = FALSE)
anno.tbl$Patient<-rownames(anno.tbl)

## Read in the clinical annotation file that lists the samples and what data
## are available for each in the columns:
data.set.cols <- colnames(anno.tbl)[grepl(pattern="File", colnames(anno.tbl)) & !grepl(pattern="Sampl", colnames(anno.tbl))]

## The mapping between the Patient/row of the clinical annotation file and
## the identifier in each data file is provided by the corresponding ID
## columns of the clinical annotation file.
## These have the same names as the data set columns, but with "SamplId"
## appended.
data.set.patient.id.cols <- unlist(lapply(data.set.cols,
                                          function(str) paste0(str, "SamplId")))


## Restrict to gene-level expression data sets
expression.data.set.flag <- grepl(pattern="geneLevelExp", data.set.cols)
data.set.cols <- data.set.cols[expression.data.set.flag]
data.set.patient.id.cols <- data.set.patient.id.cols[expression.data.set.flag]


if(use.sequestered.data) {
  ## Subset to the files we are interested in.
  all.data.sets <- unique(na.omit(as.vector(as.matrix(anno.tbl[,data.set.cols]))))
  chal_data_df <- subset(chal_data_df, name %in% all.data.sets)
  ## Download the data and put it in the base directory of the docker image,
  ## which is where it would be mounted in the challenge.
  sapply(chal_data_df$id, function(x) { synGet(x, downloadLocation=test.dir)})
}


## Read in each of the data sets and create:
## 1) A list holding each of the data sets.
## 2) A list with an entry for each data set holding a mapping from the
##    sample name in the clinical annotation column and the sample name
##    in the expression data set.
data.sets <- list()
sample.name.mappings <- list()

for(col.indx in 1:length(data.set.cols)) {
  for(data.set in na.omit(unique(anno.tbl[,data.set.cols[col.indx]]))) {
    cat(paste0("Reading in data.set ", data.set, "\n"))
    file <- paste0(test.dir, data.set)
    tbl <- as.data.frame(fread(file, header=TRUE))
    
    ## The first column is the gene identifier--make it the row
    ## name and drop it from the table
    rownames(tbl) <- tbl[,1]
    tbl <- tbl[,-1]
    data.sets[[data.set]] <- tbl
    
    ## Extract the sample name mappings from the annotation file
    patient.id.col <- data.set.patient.id.cols[col.indx]
    flag <- !is.na(anno.tbl[,data.set.cols[col.indx]]) & (anno.tbl[,data.set.cols[col.indx]] == data.set)
    map <- anno.tbl[which(flag==TRUE),c( patient.id.col)]
    map<-as.data.frame(map)
    colnames(map)<-patient.id.col
    map$Patient<-row.names(anno.tbl)[c(which(flag==TRUE))]
    map<-map[,c("Patient",patient.id.col)]
    colnames(map) <- c("patient.id", "sample.id")
    
    sample.name.mappings[[data.set]] <- map
  }
}


## Combine multiple sample columns into a single patient column by
## taking the average of the sample columns
## Columns (i.e., samples) are assumed to have names in sample.to.patient.map$from.
## They will be renamed as sample.to.patient.map$to.
combine_samples_2_patient <- function(expr, sample.to.patient.map) {
  map <- subset(sample.to.patient.map, from %in% colnames(expr))
  rownames(map) <- map$from
  expr <- expr[, colnames(expr) %in% map$from]
  map <- map[colnames(expr),]
  expr <- as.matrix(t(aggregate.Matrix(t(expr), groupings=list(map$to), fun = "mean")))
  expr
}

## Patients listed in the annotation file may have 0, 1, or >1 corresponding
## samples in the data file.  If there are multiple, take the mean across
## samples (NB: since most of these are log expression values, this corresponds
## to the geometric mean in real space).
## Also, subset both the mapping tables and the data to have the
## intersection of patients in both.
for(data.set in names(sample.name.mappings)) {
  tbl <- data.sets[[data.set]]
  map <- sample.name.mappings[[data.set]]
  
  map <- na.omit(map)
  
  ## If the mapping from patient to sample id is 1-to-many, split it into multiple
  ## rows, each of which maps the patient to a single sample id.
  map <- ldply(1:nrow(map),
               .fun = function(i) {
                 sample.id <- unlist(strsplit(as.character(map$sample.id[i]), split=";[ ]*"))
                 df <- data.frame(patient.id = map$patient.id[i], sample.id = sample.id)
                 df
               })
  
  map <- map[,c("patient.id", "sample.id")]
  colnames(map) <- c("to", "from")
  tbl <- combine_samples_2_patient(tbl, map)
  
  ## At this point, we have renamed and combined the sample ids columns to patient ids
  ## so that the map is just the identity
  map <- data.frame(patient.id = colnames(tbl), sample.id = colnames(tbl))
  
  data.sets[[data.set]] <- tbl
  sample.name.mappings[[data.set]] <- map
}

map <- do.call("rbind", sample.name.mappings)
rownames(map) <- NULL
map <- unique(map)
if(any(duplicated(map$patient.id))) {
  warning("Was not expecting any patient IDs to be duplicated\n")
}

## The sample IDs are scattered across multiple columns, put them
## in the single "sample.id" column
anno.tbl <- merge(anno.tbl, map, by.x = "Patient", by.y = "patient.id", all = FALSE)


## Subset the expression data sets to only have the genes in the genes_list
## data set (and that are in common between all data sets).
all.genes <- lapply(data.sets, rownames)

model.entrezs<-c(U_genes_list,E_genes_list,H_genes_list,M_genes_list,eu_genes_list)
model.entrezs<-unique(model.entrezs)


## Z-score each gene, separately for each data set.
## NB: scale normalizes the _columns_ of the matrix
for(data.set in names(data.sets)) {
  if (data.set=="dfci.2009_entrezID_TPM_hg19_LBR1.csv" | data.set=="m2gen_entrezID_TPM_hg19_LBR1.csv" ) {
    data.sets[[data.set]][is.na(data.sets[[data.set]])]<-0
    data.sets[[data.set]]<-Linnorm(data.sets[[data.set]])}
}
## Subset all of the data sets
for(data.set in names(data.sets)) {
  data.sets[[data.set]]<-t(data.sets[[data.set]])
  select.data<-setNames(as.data.frame(matrix(0,nrow=nrow(data.sets[[data.set]]),ncol=length(model.entrezs))),model.entrezs)
  idx<-colnames(data.sets[[data.set]])[colnames(data.sets[[data.set]])%in%colnames(select.data)]
  select.data[,idx]<-data.sets[[data.set]][,idx]
  rownames(select.data)<-rownames(data.sets[[data.set]])
  data.sets[[data.set]] <- select.data
}

expr <- do.call("rbind", data.sets)
rownames(expr)<-sub(".*.csv.","",rownames(expr))
test_MM<-expr
nfea=ncol(test_MM)

#split test data into different subset to fit in different models
test_U<-setNames(as.data.frame(matrix(0,nrow=nrow(expr),ncol=length(U_genes_list))),U_genes_list)
idx<-colnames(expr)[colnames(expr)%in%colnames(test_U)]
test_U[,idx]<-expr[,idx]
rownames(test_U)<-rownames(expr)

test_E<-setNames(as.data.frame(matrix(0,nrow=nrow(expr),ncol=length(E_genes_list))),E_genes_list)
idx<-colnames(expr)[colnames(expr)%in%colnames(test_E)]
test_E[,idx]<-expr[,idx]
rownames(test_E)<-rownames(expr)

test_H<-setNames(as.data.frame(matrix(0,nrow=nrow(expr),ncol=length(H_genes_list))),H_genes_list)
idx<-colnames(expr)[colnames(expr)%in%colnames(test_H)]
test_H[,idx]<-expr[,idx]
rownames(test_H)<-rownames(expr)

test_M<-setNames(as.data.frame(matrix(0,nrow=nrow(expr),ncol=length(M_genes_list))),M_genes_list)
idx<-colnames(expr)[colnames(expr)%in%colnames(test_M)]
test_M[,idx]<-expr[,idx]
rownames(test_M)<-rownames(expr)


test_eu<-setNames(as.data.frame(matrix(0,nrow=nrow(expr),ncol=length(eu_genes_list))),eu_genes_list)
idx<-colnames(expr)[colnames(expr)%in%colnames(test_eu)]
test_eu[,idx]<-expr[,idx]
rownames(test_eu)<-rownames(expr)

th<-0.5

##load the U_mdoel 
test_U_svm_prob<-predict(U_svmfit,test_U[,1:ncol(test_U)],probability = TRUE)
test_U_svm<-predict(U_svmfit,test_U[,1:ncol(test_U)])
test_U_rf_prob<-predict(U_rffit,test_U[,1:ncol(test_U)],type="prob")
test_U_rf<-factor(ifelse(test_U_rf_prob$X1>th,"1","-1"))
test_U_nn_prob<-predict(U_nn_fit,test_U[,1:ncol(test_U)],type="prob")
test_U_nn<-predict(U_nn_fit,test_U[,1:ncol(test_U)])
test_U_lvq<-predict(U_lvq_fit,test_U[,1:ncol(test_U)])
test_U_lvq_prob<-predict(U_lvq_fit,test_U[,1:ncol(test_U)],type="prob")
test_U_gbm_prob<-predict(U_gbm_fit,test_U[,1:ncol(test_U)],type="prob")
test_U_gbm<-factor(ifelse(test_U_gbm_prob$X1>th,"1","-1"))
test_U_glm_prob<-predict(U_glmfit,test_U[,1:ncol(test_U)],type="prob")
test_U_glm<-factor(ifelse(test_U_glm_prob$X1>th,"1","-1"))
combo_U_test<-data.frame(attributes(test_U_svm_prob)$probabilities[,2],test_U_rf_prob$X1,test_U_nn_prob$`1`,test_U_lvq_prob$X1,test_U_gbm_prob$X1,test_U_glm_prob$X1)
colnames(combo_U_test)<-c("attributes.svm_prob_train..probabilities...2.","rf_prob_train.X1","nn_prob_train..1.","lvq_prob_train.X1","gbm_prob_train.X1","glm_prob_train.X1")
test_U_prob<-predict(U_fit_stacked,combo_U_test,type="prob")
th<-0.5
test_U_stacked<-factor(ifelse(test_U_prob$X1>th,"1","-1"))


##load the E_mdoel
test_E_svm_prob<-predict(E_svmfit,test_E[,1:ncol(test_E)],probability = TRUE)
test_E_svm<-predict(E_svmfit,test_E[,1:ncol(test_E)])
test_E_rf_prob<-predict(E_rffit,test_E[,1:ncol(test_E)],type="prob")
test_E_rf<-factor(ifelse(test_E_rf_prob$X1>th,"1","-1"))
test_E_nn_prob<-predict(E_nn_fit,test_E[,1:ncol(test_E)],type="prob")
test_E_nn<-predict(E_nn_fit,test_E[,1:ncol(test_E)])
test_E_lvq<-predict(E_lvq_fit,test_E[,1:ncol(test_E)])
test_E_lvq_prob<-predict(E_lvq_fit,test_E[,1:ncol(test_E)],type="prob")
test_E_gbm_prob<-predict(E_gbm_fit,test_E[,1:ncol(test_E)],type="prob")
test_E_gbm<-factor(ifelse(test_E_gbm_prob$X1>th,"1","-1"))
test_E_glm_prob<-predict(E_glmfit,test_E[,1:ncol(test_E)],type="prob")
test_E_glm<-factor(ifelse(test_E_glm_prob$X1>th,"1","-1"))
combo_E_test<-data.frame(attributes(test_E_svm_prob)$probabilities[,2],test_E_rf_prob$X1,test_E_nn_prob$`1`,test_E_lvq_prob$X1,test_E_gbm_prob$X1,test_E_glm_prob$X1)
colnames(combo_E_test)<-c("attributes.svm_prob_train..probabilities...2.","rf_prob_train.X1","nn_prob_train..1.","lvq_prob_train.X1","gbm_prob_train.X1","glm_prob_train.X1")
test_E_prob<-predict(E_fit_stacked,combo_E_test,type="prob")
th<-0.5
test_E_stacked<-factor(ifelse(test_E_prob$X1>th,"1","-1"))

##load the H_mdoel
test_H_svm_prob<-predict(H_svmfit,test_H[,1:ncol(test_H)],probability = TRUE)
test_H_svm<-predict(H_svmfit,test_H[,1:ncol(test_H)])
test_H_rf_prob<-predict(H_rffit,test_H[,1:ncol(test_H)],type="prob")
test_H_rf<-factor(ifelse(test_H_rf_prob$X1>th,"1","-1"))
test_H_nn_prob<-predict(H_nn_fit,test_H[,1:ncol(test_H)],type="prob")
test_H_nn<-predict(H_nn_fit,test_H[,1:ncol(test_H)])
test_H_lvq<-predict(H_lvq_fit,test_H[,1:ncol(test_H)])
test_H_lvq_prob<-predict(H_lvq_fit,test_H[,1:ncol(test_H)],type="prob")
test_H_gbm_prob<-predict(H_gbm_fit,test_H[,1:ncol(test_H)],type="prob")
test_H_gbm<-factor(ifelse(test_H_gbm_prob$X1>th,"1","-1"))
test_H_glm_prob<-predict(H_glmfit,test_H[,1:ncol(test_H)],type="prob")
test_H_glm<-factor(ifelse(test_H_glm_prob$X1>th,"1","-1"))
combo_H_test<-data.frame(attributes(test_H_svm_prob)$probabilities[,2],test_H_rf_prob$X1,test_H_nn_prob$`1`,test_H_lvq_prob$X1,test_H_gbm_prob$X1,test_H_glm_prob$X1)
colnames(combo_H_test)<-c("attributes.svm_prob_train..probabilities...2.","rf_prob_train.X1","nn_prob_train..1.","lvq_prob_train.X1","gbm_prob_train.X1","glm_prob_train.X1")
test_H_prob<-predict(H_fit_stacked,combo_H_test,type="prob")
th<-0.5
test_H_stacked<-factor(ifelse(test_H_prob$X1>th,"1","-1"))



##load the M_mdoel
test_M_svm<-predict(M_svmfit,test_M[,1:ncol(test_M)])
test_M_svm_prob<-predict(M_svmfit,test_M[,1:ncol(test_M)],probability = TRUE)
test_M_rf_prob<-predict(M_rffit,test_M[,1:ncol(test_M)],type="prob")
test_M_rf<-factor(ifelse(test_M_rf_prob$X1>th,"1","-1"))
test_M_nn_prob<-predict(M_nn_fit,test_M[,1:ncol(test_M)],type="prob")
test_M_nn<-predict(M_nn_fit,test_M[,1:ncol(test_M)])
test_M_lvq<-predict(M_lvq_fit,test_M[,1:ncol(test_M)])
test_M_lvq_prob<-predict(M_lvq_fit,test_M[,1:ncol(test_M)],type="prob")
test_M_gbm_prob<-predict(M_gbm_fit,test_M[,1:ncol(test_M)],type="prob")
test_M_gbm<-factor(ifelse(test_M_gbm_prob$X1>th,"1","-1"))
test_M_glm_prob<-predict(M_glmfit,test_M[,1:ncol(test_M)],type="prob")
test_M_glm<-factor(ifelse(test_M_glm_prob$X1>th,"1","-1"))
th<-0.5
combo_M_test<-data.frame(attributes(test_M_svm_prob)$probabilities[,2],test_M_rf_prob$X1,test_M_nn_prob$`1`,test_M_lvq_prob$X1,test_M_gbm_prob$X1,test_M_glm_prob$X1)
colnames(combo_M_test)<-c("attributes.svm_prob_train..probabilities...2.","rf_prob_train.X1","nn_prob_train..1.","lvq_prob_train.X1","gbm_prob_train.X1","glm_prob_train.X1")
test_M_prob<-predict(M_fit_stacked,combo_M_test,type="prob")
th<-0.5
test_M_stacked<-factor(ifelse(test_M_prob$X1>th,"1","-1"))


##load the eu_mdoel 
test_eu_svm<-predict(eu_svmfit,test_eu[,1:ncol(test_eu)])
test_eu_svm_prob<-predict(eu_svmfit,test_eu[,1:ncol(test_eu)],probability = TRUE)
test_eu_rf_prob<-predict(eu_rffit,test_eu[,1:ncol(test_eu)],type="prob")
test_eu_rf<-factor(ifelse(test_eu_rf_prob$X1>th,"1","-1"))
test_eu_nn_prob<-predict(eu_nn_fit,test_eu[,1:ncol(test_eu)],type="prob")
test_eu_nn<-predict(eu_nn_fit,test_eu[,1:ncol(test_eu)])
test_eu_lvq<-predict(eu_lvq_fit,test_eu[,1:ncol(test_eu)])
test_eu_lvq_prob<-predict(eu_lvq_fit,test_eu[,1:ncol(test_eu)],type="prob")
test_eu_gbm_prob<-predict(eu_gbm_fit,test_eu[,1:ncol(test_eu)],type="prob")
test_eu_gbm<-factor(ifelse(test_eu_gbm_prob$X1>th,"1","-1"))
test_eu_glm_prob<-predict(eu_glmfit,test_eu[,1:ncol(test_eu)],type="prob")
test_eu_glm<-factor(ifelse(test_eu_glm_prob$X1>th,"1","-1"))
th<-0.5
combo_eu_test<-data.frame(attributes(test_eu_svm_prob)$probabilities[,2],test_eu_rf_prob$X1,test_eu_nn_prob$`1`,test_eu_lvq_prob$X1,test_eu_gbm_prob$X1,test_eu_glm_prob$X1)
colnames(combo_eu_test)<-c("attributes.svm_prob_train..probabilities...2.","rf_prob_train.X1","nn_prob_train..1.","lvq_prob_train.X1","gbm_prob_train.X1","glm_prob_train.X1")
test_eu_prob<-predict(eu_fit_stacked,combo_eu_test,type="prob")
th<-0.5
test_eu_stacked<-factor(ifelse(test_eu_prob$X1>th,"1","-1"))


##global ensemble_test data
U_models_test<-cbind(combo_U_test[,1:6],test_U_prob$X1)
colnames(U_models_test)<-gsub("_prob_","_U_",colnames(U_models_test))
colnames(U_models_test)[7]<-c("train_U_prob$X1")
E_models_test<-cbind(combo_E_test[,1:6],test_E_prob$X1)
colnames(E_models_test)<-gsub("_prob_","_E_",colnames(E_models_test))
colnames(E_models_test)[7]<-c("train_E_prob$X1")
H_models_test<-cbind(combo_H_test[,1:6],test_H_prob$X1)
colnames(H_models_test)<-gsub("_prob_","_H_",colnames(H_models_test))
colnames(H_models_test)[7]<-c("train_H_prob$X1")
M_models_test<-cbind(combo_M_test[,1:6],test_M_prob$X1)
colnames(M_models_test)<-gsub("_prob_","_M_",colnames(M_models_test))
colnames(M_models_test)[7]<-c("train_M_prob$X1")
eu_models_test<-cbind(combo_eu_test[,1:6],test_eu_prob$X1)
colnames(eu_models_test)<-gsub("_prob_","_eu_",colnames(eu_models_test))
colnames(eu_models_test)[7]<-c("train_eu_prob$X1")

ensemble_test<-cbind(U_models_test,E_models_test,H_models_test,M_models_test,eu_models_test)
ensemble_test$label<-combo_U_test$train_label
colnames(ensemble_test)<-ensemble_rf$coefnames
colnames(ensemble_test)<-c("attributes.svm_U_train..probabilities...2.","rf_U_train.X1",                              
                           "nn_U_train..1.","lvq_U_train.X1",                             
                           "gbm_U_train.X1","glm_U_train.X1",                             
                           "train_U_prob$X1","attributes.svm_E_train..probabilities...2.", 
                           "rf_E_train.X1","nn_E_train..1.",                             
                           "lvq_E_train.X1","gbm_E_train.X1",                             
                           "glm_E_train.X1","train_E_prob$X1",                            
                           "attributes.svm_H_train..probabilities...2.",  "rf_H_train.X1",                              
                           "nn_H_train..1."  ,                            "lvq_H_train.X1",                             
                           "gbm_H_train.X1"  ,                            "glm_H_train.X1" ,                            
                           "train_H_prob$X1"  ,                           "attributes.svm_M_train..probabilities...2." ,
                           "rf_M_train.X1"  ,                             "nn_M_train..1.",                             
                           "lvq_M_train.X1" ,                             "gbm_M_train.X1" ,                            
                           "glm_M_train.X1" ,                             "train_M_prob$X1"  ,                          
                           "attributes.svm_eu_train..probabilities...2." ,"rf_eu_train.X1" ,                            
                           "nn_eu_train..1." ,                            "lvq_eu_train.X1"   ,                         
                           "gbm_eu_train.X1" ,                            "glm_eu_train.X1" ,                           
                           "train_eu_prob$X1"  )




ensemble_rf_prob<-predict(ensemble_rf,ensemble_test[,1:35],type="prob")  ###
th<-0.5
ensemble_rf_pred<-factor(ifelse(ensemble_rf_prob$X1>th,"1","-1"))

###generate prediction score
predictionscore<-array()
PS_U<-array()
PS_E<-array()
PS_H<-array()
PS_M<-array()
for (i in 1:length(ensemble_rf_pred)){
  PS_U[i]<-sum(attr(test_U_svm_prob,"probabilities")[i,2],test_U_gbm_prob$X1[i],test_U_glm_prob$X1[i],test_U_nn_prob[i,"1"],test_U_rf_prob[i,"1"],test_U_prob$X1[i],test_U_lvq_prob$X1[i])
  PS_H[i]<-sum(attr(test_H_svm_prob,"probabilities")[i,2],test_H_gbm_prob$X1[i],test_H_glm_prob$X1[i],test_H_nn_prob[i,"1"],test_H_rf_prob[i,"1"],test_H_prob$X1[i],test_H_lvq_prob$X1[i])
  PS_E[i]<-sum(attr(test_E_svm_prob,"probabilities")[i,2],test_E_gbm_prob$X1[i],test_E_glm_prob$X1[i],test_E_nn_prob[i,"1"],test_E_rf_prob[i,"1"],test_E_prob$X1[i],test_E_lvq_prob$X1[i])
  PS_M[i]<-sum(attr(test_M_svm_prob,"probabilities")[i,2],test_M_gbm_prob$X1[i],test_M_glm_prob$X1[i],test_M_nn_prob[i,"1"],test_M_rf_prob[i,"1"],test_M_prob$X1[i],test_M_lvq_prob$X1[i])
  if (ensemble_rf_pred[i]=="1"){
    predictionscore[i]<-max(c(PS_U[i],PS_H[i],PS_E[i]),PS_M[i])
  } else {
    predictionscore[i]<-min(c(PS_U[i],PS_H[i],PS_E[i]),PS_M[i])
  }
}

predictionscore<-round(predictionscore,5)
highriskflag<-as.character(ensemble_rf_pred)
for (i in 1 : length(highriskflag)) {if (highriskflag[i]=="-1") {highriskflag[i]="False"}else {highriskflag[i]="True"} }
patient<-rownames(ensemble_test)
S_P<-data.frame(anno.tbl$Study,anno.tbl$Patient)
colnames(S_P)<-c("study","patient")


## Create the output directory
dir.create(output.dir)
## Output the results with columns ID, raw.score, and high.risk.
output.file <- paste0(output.dir, "/", "predictions.tsv")
output1<-data.frame(patient,predictionscore,highriskflag)
output<-merge(S_P,output1,by="patient",all.x=TRUE)
predictions.tbl<-output[,c("study","patient","predictionscore","highriskflag")]
predictions.tbl$patient<-gsub('\\.','_',predictions.tbl$patient)
write.table(file=output.file, predictions.tbl,row.names = FALSE,col.names=TRUE, quote=FALSE, sep="\t")

cat("Successfully wrote predictions.\n")

