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
  chal_data_df_anno <- subset(chal_data_df, grepl(name, pattern="sc3_Validation_ClinAnnotations.csv"))
  sapply(chal_data_df_anno$id, function(x) { synGet(x, downloadLocation=test.dir)})
}


## Read in the annotation file
training.validation.file <- paste0(test.dir, "sc3_Validation_ClinAnnotations.csv")
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

model.entrezs<-c(M_genes_eu_list)
model.entrezs<-unique(model.entrezs)


## Normalize RNA_seq data
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
for (i in 1: ncol(expr)){
  expr[,i][is.na(expr[,i])] <- median(as.numeric(expr[,i]),na.rm=T)
}
rownames(expr)<-sub(".*.csv.","",rownames(expr))

test_MM<-expr
th<-0.5
test_M<-test_MM

risk_flag_M<-setNames(as.data.frame(matrix(0,ncol=2,nrow = nrow(test_M))),c("risk_flag","predictionscore"))
rownames(risk_flag_M)<-rownames(test_M)

Age_ISS<-as.data.frame(anno.tbl[,c("D_Age","D_ISS")])
for (i in 1: ncol(Age_ISS)){
  Age_ISS[,i][is.na(Age_ISS[,i])] <- median(as.numeric(Age_ISS[,i]),na.rm=T)
}
rownames(Age_ISS)<-anno.tbl$Patient
Age_ISS$D_Age<-Age_ISS$D_Age/10


#M_model
test_M_EU<-setNames(as.data.frame(matrix(0,nrow=nrow(test_M),ncol=length(M_genes_eu_list))),M_genes_eu_list)
rownames(test_M_EU)<-rownames(test_M)
idx<-colnames(test_M)[colnames(test_M)%in%colnames(test_M_EU)]
test_M_EU[,idx]<-test_M[,idx]
test_M<-test_M_EU
##load the M_mdoel
test_M_svm_prob<-predict(M_svmfit,test_M[,1:ncol(test_M)],probability = TRUE)
test_M_svm<-predict(M_svmfit,test_M[,1:ncol(test_M)])
test_M_rf_prob<-predict(M_rffit,test_M[,1:ncol(test_M)],type="prob")
test_M_rf<-factor(ifelse(test_M_rf_prob$X1>th,"1","-1"))
test_M_nn_prob<-predict(M_nn_fit,test_M[,1:ncol(test_M)],type="prob")
test_M_nn<-predict(M_nn_fit,test_M[,1:ncol(test_M)])
test_M_lvq_prob<-predict(M_lvq_fit,test_M[,1:ncol(test_M)],type="prob")
test_M_lvq<-factor(ifelse(test_M_lvq_prob$X1>th,"1","-1"))
test_M_gbm_prob<-predict(M_gbm_fit,test_M[,1:ncol(test_M)],type="prob")
test_M_gbm<-factor(ifelse(test_M_gbm_prob$X1>th,"1","-1"))
test_M_glm_prob<-predict(M_glmfit,test_M[,1:ncol(test_M)],type="prob")
test_M_glm<-factor(ifelse(test_M_glm_prob$X1>th,"1","-1"))
test_M_knn_prob<-predict(M_knnfit,Age_ISS[,1:2],type='prob')
test_M_knn<-predict(M_knnfit,Age_ISS[,1:2])
combo_M_test<-data.frame(test_M_svm,test_M_rf,test_M_nn,test_M_lvq,test_M_gbm,test_M_glm,test_M_knn)
colnames(combo_M_test)<-c( "svm_pred_train" ,"rf_pred_train",  "nn_pred_train" , "lvq_pred_train", "gbm_pred_train" ,"glm_pred_train","knn_pred_train")
test_M_prob<-predict(M_fit_stacked,combo_M_test,type="prob")
test_M_stacked<-factor(ifelse(test_M_prob$X1>th,"1","-1"))
for (i in 1: length(test_M_stacked)) {if (test_M_stacked[i]=="1"){test_M_stacked[i]=1} else {test_M_stacked[i]=-1}}
risk_flag_M[,1]<-test_M_stacked
for (i in 1: nrow(risk_flag_M)) {
  risk_flag_M$predictionscore[i]<-sum(attr(test_M_svm_prob,"probabilities")[i,2],test_M_gbm_prob$X1[i],test_M_glm_prob$X1[i],test_M_nn_prob[i,"1"],test_M_rf_prob$X1[i],test_M_lvq_prob$X1[i],test_M_prob$X1[i],test_M_knn_prob$`1`[i])
}


risk_flag<-risk_flag_M
colnames(risk_flag)<-c("highriskflag","predictionscore")
risk_flag<-risk_flag[,c('predictionscore',"highriskflag")]
risk_flag$highriskflag<-as.character(risk_flag$highriskflag)
for (i in 1 : nrow(risk_flag)) {if (risk_flag$highriskflag[i]=="-1") {risk_flag$highriskflag[i]="False"}else {risk_flag$highriskflag[i]="True"} }
patient<-rownames(risk_flag)
S_P<-data.frame(anno.tbl$Study,anno.tbl$Patient)
colnames(S_P)<-c("study","patient")


## Output the results with columns ID, raw.score, and high.risk.
output.file <- paste0(output.dir, "/", "predictions.tsv")
output1<-data.frame(patient,risk_flag)
output<-merge(S_P,output1,by="patient",all.x=TRUE)
predictions.tbl<-output[,c("study","patient","predictionscore","highriskflag")]
predictions.tbl$patient<-gsub("\\.","_",predictions.tbl$patient)
write.table(file=output.file, predictions.tbl,row.names = FALSE,col.names=TRUE, quote=FALSE, sep="\t")

cat("Successfully wrote predictions.\n")


