
#----
#  title:'Multiple Myeloma DREAM Challenge'
#  author:'Haibin Guan'
#  date:'16 July 2017'
#----

#There are four parts to my script as follows:
#    differential expression
#    clustering
#    feature engineering
#    prediction
    
#load packages
source("http://bioconductor.org/biocLite.R")
biocLite("data.table")
biocLite("limma")
biocLite("Biobase")
biocLite("genefilter")

library('ggplot2') #visualization
library('ggthemes') #visualization
library('scales') #visualization
library(data.table)
library(limma)
library(Biobase)
library(genefilter)
library(data.table)
# set working directory
setwd("~/Documents/IT Course Notes/MM/dataset/DREAM_UAMS_data")

#read files in r
Gene_expression=read.csv(file="GSE24080UAMSentrezIDlevel.csv",header=T,sep=",")
clinical<-read.csv("globalClinTraining.csv",header=T,sep=",",row.names=2)
clinical<-subset(clinical,Study=='GSE24080UAMS')
clinical<-subset(clinical[,2:8])
pData<-AnnotatedDataFrame(clinical)

#convert file to a ExpressionSet
ex<-new("ExpressionSet",exprs=as.matrix(Gene_expression),phenoData=pData)

######
# > ex
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 20514 features, 559 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM592391 GSM592392 ... GSM592661 (559 total)
# varLabels: D_Age D_Gender ... D_ISS (7 total)
# varMetadata: labelDescription
# featureData: none
# experimentData: use 'experimentData(object)'
# Annotation:

######


library(genefilter)
ex_sd<-rowSds(exprs(ex))
top5000<-names(sort(ex_sd,decreasing=T))[1:5000]
ex_var<-ex[top5000,]
dist_cor <- function(x) {
       as.dist(1 - cor(t(x), method = "pearson"))
   }
 clus_wd2 <- function(x) {
      hclust(x, method = "ward.D2")
   }
 library(RColorBrewer)
 redblackgreen <- colorRampPalette(c("green", "black", "red"))(n = 100)

 class_labels <- ifelse(ex_var$D_PFS_FLAG == "1", "grey80", "grey20")

   library(gplots)
 heatmap.2(exprs(ex_var), 
                        # clustering
                          distfun = dist_cor, 
                        hclust = clus_wd2,
                        # scaling (genes are in rows)
                          scale = "row",
                        # color
                          col = redblackgreen, 
                        # labels
                          labRow = "", 
                        # tweaking
                          trace = "none",
                        density.info = "none")


 ##code modified based on 'deepnet' package
 ##Training a Deep neural network with weights initialized by Stacked AutoEncoder
 sae.dnn.train <- function(x,y,hidden=c(500,100,50),
                           activationfun="sigm",  #activationfun activation function of hidden unit. "sigm" for logistic function
                           learningrate=0.5,   #learningrate learning rate for gradient descent.
                           momentum=0.3,
                           learningrate_scale=0.1, #learningrate_scale  learning rate will be mutiplied by this scale after every iteration.
                           output="sigm",
                           sae_output="linear",
                           numepochs=3,batchsize=100, #number of iteration for samples; size of mini-batch
                           hidden_dropout=0,visible_dropout=0){
   
   message("begin to train sae ......")
   sae <- sae.train(x,hidden=hidden,
                    activationfun=activationfun,
                    output=sae_output,
                    numepochs=numepochs,batchsize=batchsize,
                    learningrate=learningrate,learningrate_scale=learningrate_scale,
                    momentum=momentum,
                    hidden_dropout=hidden_dropout,visible_dropout=visible_dropout)
   message("sae has been trained.")
   initW <- list()
   initB <- list()
   for(i in 1:(length(sae$size) - 1)){
     initW[[i]] <- sae$encoder[[i]]$W[[1]]   
     initB[[i]] <- sae$encoder[[i]]$B[[1]]  
   }
   
   
   #random init weights between last hidden layer and output layer
   last_hidden <- sae$size[length(sae$size)]
   initW[[length(sae$size)]] <- matrix(runif(output_dim*last_hidden,min=-0.1,max=0.1), c(output_dim,last_hidden))
   initB[[length(sae$size)]] <- runif(output_dim,min=-0.1,max=0.1)
   message("begin to train deep nn ......")
   dnn <- nn.train(x,y,initW=initW,initB=initB,hidden=hidden,
                   activationfun=activationfun,
                   learningrate=learningrate,
                   momentum=momentum,
                   learningrate_scale=learningrate_scale,
                   output=output,
                   numepochs=numepochs,batchsize=batchsize,
                   hidden_dropout=hidden_dropout,visible_dropout=visible_dropout)
   message("deep nn has been trained.")
   dnn
 }
 
 
 
 
 sae.train <- function(x,hidden=c(10),
                       activationfun="sigm",
                       learningrate=0.8,
                       momentum=0.5,
                       learningrate_scale=1,
                       output="sigm",
                       numepochs=3,batchsize=100,
                       hidden_dropout=0,visible_dropout=0.2                      
 ){
   input_dim <- ncol(x)
   size <- c(input_dim, hidden)
   sae <- list(
     input_dim = input_dim,
     hidden = hidden,
     size = size
   )
   train_x <- x
   message(sprintf("training layer 1 autoencoder ..."))
   sae$encoder[[1]] <-  nn.train(train_x,train_x,hidden=c(hidden[1]),
                                 activationfun=activationfun,
                                 learningrate=learningrate,
                                 momentum=momentum,
                                 learningrate_scale=learningrate_scale,
                                 output=output,
                                 numepochs=numepochs,batchsize=batchsize,
                                 hidden_dropout=hidden_dropout,visible_dropout=visible_dropout)
   
   if(length(sae$size) > 2){
     for(i in 2:(length(sae$size) - 1)){
       pre <- t( sae$encoder[[i-1]]$W[[1]] %*% t(train_x) + sae$encoder[[i-1]]$B[[i-1]] )
       if(sae$encoder[[i-1]]$activationfun == "sigm"){
         post <- sigm( pre )
       }else if(sae$encoder[[i-1]]$activationfun == "tanh"){
         post <- tanh(pre)
       }else{
         stop("unsupport activation function 'nn$activationfun'");
       }  
       train_x <- post
       message(sprintf("training layer %d autoencoder ...",i))
       sae$encoder[[i]] <- nn.train(train_x,train_x,hidden=c(hidden[i]),
                                    activationfun=activationfun,
                                    learningrate=learningrate,
                                    momentum=momentum,
                                    learningrate_scale=learningrate_scale,
                                    output=output,
                                    numepochs=numepochs,batchsize=batchsize,
                                    hidden_dropout=hidden_dropout,visible_dropout=visible_dropout)
     }
   }
   sae
 }
 
 
 
 
 nn.train <- function(x,y,initW=NULL,initB=NULL,hidden=c(10),
                      activationfun="sigm",
                      learningrate=0.8,
                      momentum=0.5,
                      learningrate_scale=1,
                      output="sigm",
                      numepochs=3,batchsize=100,
                      hidden_dropout=0,visible_dropout=0) {
   
   size <- c(input_dim, hidden, output_dim)
   vW <- list() 
   vB <- list()
   if(is.null(initW) || is.null(initB)){
     W <- list()
     B <- list()
     #random init weights and bias between layers							 
     for( i in 2:length(size) ){
       W[[i-1]] <- matrix(runif(size[i]*size[i-1],min=-0.1,max=0.1), c(size[i],size[i-1]));
       B[[i-1]] <- runif(size[i],min=--0.1,max=0.1);
       vW[[i-1]] <- matrix(rep(0,size[i]*size[i-1]),c(size[i],size[i-1]))
       vB[[i-1]] <- rep(0,size[i])
     }
   }else{
     W <- initW
     B <- initB
     for( i in 2:length(size) ){
       vW[[i-1]] <- matrix(rep(0,size[i]*size[i-1]),c(size[i],size[i-1]))
       vB[[i-1]] <- rep(0,size[i])
       if(nrow(W[[i-1]]) != size[i] || ncol(W[[i-1]]) != size[i-1] ){
         stop("init W size is not eq to network size!")  
       }    
       if(length(B[[i-1]]) != size[i]){
         stop("init B size is not eq to network size!")  
       }
     }
   }
   
   nn <- list(
     input_dim = input_dim,
     output_dim = output_dim,
     hidden = hidden,
     size = size,
     activationfun = activationfun,
     learningrate = learningrate,
     momentum = momentum,
     learningrate_scale = learningrate_scale,
     hidden_dropout=hidden_dropout,visible_dropout=visible_dropout,
     output = output,
     W = W,
     vW = vW,
     B = B,
     vB = vB
   )
   
   m <- nrow(x);
   numbatches <- m / batchsize;
   s <- 0
   for(i in 1:numepochs){
     randperm <- sample(1:m,m)
     if(numbatches >= 1){
       for(l in 1 : numbatches){
         s <- s + 1
         batch_x <- x[randperm[((l-1)*batchsize+1):(l*batchsize)], ] 
         if(is.vector(y)){
           batch_y <- y[randperm[((l-1)*batchsize+1):(l*batchsize)]] 
         }else if(is.matrix(y)){
           batch_y <- y[randperm[((l-1)*batchsize+1):(l*batchsize)], ] 
         }
         
         nn <- nn.ff(nn,batch_x,batch_y,s)
         nn <- nn.bp(nn)						
       }
     }
     #last fraction of sample
     if(numbatches > as.integer(numbatches)){
       batch_x <- x[randperm[(as.integer(numbatches)*batchsize):m], ]
       if(is.vector(y)){
         batch_y <- y[randperm[(as.integer(numbatches)*batchsize):m]]
       }else if(is.matrix(y)){
         batch_y <- y[randperm[(as.integer(numbatches)*batchsize):m], ]
       }
       s <- s + 1
       nn <- nn.ff(nn,batch_x,batch_y,s)
       nn <- nn.bp(nn)	
     }
     
     nn$learningrate <- nn$learningrate * nn$learningrate_scale;
   }
   
   nn
 }
 
 nn.ff <- function(nn,batch_x,batch_y,s){
   m <- nrow(batch_x)
   #do input dropout
   if(nn$visible_dropout > 0){
     nn$dropout_mask[[1]] <- dropout.mask(ncol(batch_x),nn$visible_dropout)
     batch_x <- t( t(batch_x) * nn$dropout_mask[[1]] )
   }
   nn$post[[1]] <- batch_x
   for(i in 2:(length(nn$size) - 1)){
     nn$pre[[i]] <- t( nn$W[[i-1]] %*% t(nn$post[[(i-1)]])  + nn$B[[i-1]] )
     if(nn$activationfun == "sigm"){
       nn$post[[i]] <- sigm(nn$pre[[i]])
     }else if(nn$activationfun == "tanh"){
       nn$post[[i]] <- tanh(nn$pre[[i]])
     }else{
       stop("unsupport activation function!");
     }	
     if(nn$hidden_dropout > 0){
       nn$dropout_mask[[i]] <- dropout.mask(ncol(nn$post[[i]]),nn$hidden_dropout)
       nn$post[[i]] <- t( t(nn$post[[i]]) * nn$dropout_mask[[i]] )
     }
   }
   #output layer
   i <- length(nn$size)
   nn$pre[[i]] <- t( nn$W[[i-1]] %*% t(nn$post[[(i-1)]])  + nn$B[[i-1]] )
   if(nn$output == "sigm"){
     nn$post[[i]] <- sigm(nn$pre[[i]])
     nn$e <- batch_y - nn$post[[i]]
     nn$L[ s ] <- 0.5*sum(nn$e^2)/m
   }else if(nn$output == "linear"){
     nn$post[[i]] <- nn$pre[[i]]
     nn$e <- batch_y - nn$post[[i]]
     nn$L[ s ] <- 0.5*sum(nn$e^2)/m
   }else if(nn$output == "softmax"){
     nn$post[[i]] <- exp(nn$pre[[i]])
     nn$post[[i]] <- nn$post[[i]] / rowSums(nn$post[[i]]) 
     nn$e <- batch_y - nn$post[[i]]
     nn$L[ s ] <- -sum(batch_y * log(nn$post[[i]]))/m
   }else{
     stop("unsupport output function!");
   }	
   if(s %% 10000 == 0){
     message(sprintf("####loss on step %d is : %f",s,nn$L[ s ]))
   }
   
   nn
 }
 
 
 nn.bp <- function(nn){
   n <- length(nn$size)
   d <- list()
   if(nn$output == "sigm"){
     d[[n]] <- -nn$e * (nn$post[[n]] * (1 - nn$post[[n]]))
   }else if(nn$output == "linear" || nn$output == "softmax"){
     d[[n]] <- -nn$e
   }
   
   for( i in (n-1):2 ){
     if(nn$activationfun  == "sigm"){
       d_act <- nn$post[[i]] * (1-nn$post[[i]])
     }else if(nn$activationfun  == "tanh" ){
       d_act <- 1.7159 * 2/3 * (1 - 1/(1.7159)^2 * nn$post[[i]]^2)
     }
     d[[i]] <- (d[[i+1]] %*% nn$W[[i]]) * d_act
     if(nn$hidden_dropout > 0){
       d[[i]] <- t( t(d[[i]]) * nn$dropout_mask[[i]] )
     }
   }
   
   for( i in 1:(n-1) ){
     dw <- t(d[[i+1]]) %*% nn$post[[i]] / nrow(d[[i+1]])
     dw <- dw * nn$learningrate
     if(nn$momentum > 0){
       nn$vW[[i]] <- nn$momentum * nn$vW[[i]] + dw
       dw <- nn$vW[[i]]
     }
     nn$W[[i]] <- nn$W[[i]] - dw
     
     db <- colMeans(d[[i+1]])
     db <- db * nn$learningrate
     if(nn$momentum > 0){
       nn$vB[[i]] <- nn$momentum * nn$vB[[i]] + db
       db <- nn$vB[[i]]
     }
     nn$B[[i]] <- nn$B[[i]] - db
   }
   nn
 }
 
 dropout.mask <- function(size,fraction){
   mask <- runif(size,min=0,max=1)
   mask[mask <= fraction] <- 0
   mask[mask > fraction] <- 1
   mask
 }
 
 
 sigm <- function(x){
   1/(1+exp(-x))
 }
 
 dnn<-dbn.dnn.train(train_x, train_y, hidden=c(500,100,50), activationfun="sigm", learningrate=0.7, momentum=0.6, learningrate_scale=1, output="sigm",numepochs=3, batchsize=50, hidden_dropout=0, visible_dropout=0, cd=1)
 
 
 
 
 
 ##SAM
 train_x<-t(train_x)
 data=list(x=train_x,y=train_y,geneid=rownames(train_x),genenames=as.character(rownames(train_x),logged2=TRUE))
 

 
 
 
 ##Deep Learning in H2o
 dl_fit1 <- h2o.deeplearning(x = x,
                             y = y,
                             training_frame = train,
                             model_id = "dl_fit1",
                             hidden = c(20,20),
                             seed = 1)
 
 dl_fit2 <- h2o.deeplearning(x = x,
                             y = y,
                             training_frame = train,
                             model_id = "dl_fit2",
                             epochs = 50,
                             hidden = c(20,20),
                             stopping_rounds = 0,  # disable early stopping
                             seed = 1)
 
 dl_fit3 <- h2o.deeplearning(x = x,
                             y = y,
                             training_frame = train,
                             model_id = "dl_fit3",
                             epochs = 50,
                             hidden = c(20,20),
                             nfolds = 3,                            #used for early stopping
                             score_interval = 1,                    #used for early stopping
                             stopping_rounds = 5,                   #used for early stopping
                             stopping_metric = "misclassification", #used for early stopping
                             stopping_tolerance = 1e-3,             #used for early stopping
                             seed = 1)
 dl_perf1 <- h2o.performance(model = dl_fit1, newdata = test)
 dl_perf2 <- h2o.performance(model = dl_fit2, newdata = test)
 dl_perf3 <- h2o.performance(model = dl_fit3, newdata = test)
 
 # Retreive test set MSE
 h2o.mse(dl_perf1)
 
 h2o.scoreHistory(dl_fit3)
 
 h2o.confusionMatrix(dl_fit3)

 # “plot a model”, which will graph the performance of some metric over the training proces 
 plot(dl_fit3, 
      timestep = "epochs", 
      metric = "classification_error")
 
 # Get the CV models from the `dl_fit3` object
 cv_models <- sapply(dl_fit3@model$cross_validation_models, 
                     function(i) h2o.getModel(i$name))
 
 # Plot the scoring history over time
 plot(cv_models[[1]], 
      timestep = "epochs", 
      metric = "classification_error")
 
 
 ##Differential Expression Analysis: SAM 
 library(impute)
 library(matrixStats)
 library(samr)
 #create the samr data file sam.data as a list
 sam.data<-list(x=t(newdat[,1:20514]),y=t(clinical$HR_FLAG),geneid=rownames(Gene_expression),logged2=TRUE)
 samr.obj<-samr(sam.data,resp.type="Two class unpaired",nperm=616)
 names(samr.obj)
 delta.table<-samr.compute.delta.table(samr.obj,min.foldchange=1.5)
 delta.table
 delta<-0.05
 samr.plot(samr.obj,delta) ##failed Error in plot.new(): figure margins too large
 siggenes.table<-samr.compute.siggenes.table(samr.obj,delta,sam.data,delta.table )
 siggenes.table
 myupgenes<-siggenes.table$genes.up
 tail(myupgenes)
 mydowngenes<-siggenes.table$genes.lo
 mydowngenes
 ### p_value are 0s ???
 
 
 ##differential expression with genefilter package
 library(Biobase)
 library(genefilter)
 phdata<-as.data.frame(newdat$`data$Y`)
 row.names(phdata)<-row.names(newdat)
 pData<-AnnotatedDataFrame(phdata)
 
 #convert file to a ExpressionSet
 ex<-new("ExpressionSet",exprs=as.matrix(t(newdat[,1:20514])),phenoData=pData)
 
 ######
 # > ex
 # ExpressionSet (storageMode: lockedEnvironment)
 # assayData: 20514 features, 559 samples 
 # element names: exprs 
 # protocolData: none
 # phenoData
 # sampleNames: GSM592391 GSM592392 ... GSM592661 (559 total)
 # varLabels: D_Age D_Gender ... D_ISS (7 total)
 # varMetadata: labelDescription
 # featureData: none
 # experimentData: use 'experimentData(object)'
 # Annotation:
 
 ######
 library(genefilter)
 ex_sd<-rowSds(exprs(ex))
 top5000<-names(sort(ex_sd,decreasing=T))[1:5000]
 ex_var<-ex[top5000,]
 genefilter_data<-as.data.frame(t(exprs(ex_var)))
 genefilter_data$label<-pData(ex_var)
 write.csv(genefilter_data,file="genefilter_data.csv")
 
 
 
 ##visualization of features
 library(GGally)
 ggpairs(deep.fea,axisLabels = "none",mapping = ggplot2::aes(color=label,alpha=0.5),upper=list(continuous="cor"),lower=list(continuous="points"))
 
 
 # weights for different layers
 #W1<-as.data.frame(h2o.weights(MM.ae,1))
 #W2<-as.data.frame(h2o.weights(MM.ae,2))
 #W3<-as.data.frame(h2o.weights(MM.ae,3))
 #W4<-as.data.frame(h2o.weights(MM.ae,4))
 
 
 
 
 
 
 
