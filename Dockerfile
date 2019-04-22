FROM r-base:3.3.3
RUN apt-get update
RUN apt-get -y install r-base
RUN apt-get -y install aptitude libcurl4-openssl-dev libxml2-dev libxml2-dev 

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite("plyr")'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite("Matrix.utils")'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite("Linnorm")'
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('data.table')"
RUN Rscript -e "install.packages('FCNN4R')"
RUN Rscript -e "install.packages('caret')"
RUN Rscript -e "install.packages('ranger')"
RUN Rscript -e "install.packages('e1071')"
RUN Rscript -e "install.packages('randomForest')"
RUN Rscript -e "install.packages('xgboost')"
RUN Rscript -e "install.packages('glmnet')"
RUN Rscript -e "install.packages('gbm')"
RUN Rscript -e "install.packages('class')"
RUN Rscript -e "install.packages('import')"
RUN Rscript -e "install.packages('matrixStats')"



COPY score_sc2.sh /score_sc2.sh
COPY run-mm-sc2.R /run-mm-sc2.R
RUN chmod 770 /score_sc2.sh
COPY model-state-metadata.Rd /model-state-metadata.Rd