# When using, please cite: Olave, M., Avila, L.J., Sites, J.W. Jr and Morando, M (2017). Detecting Hybridization by Likelihood Calculation of Gene Tree Extra Lineages given Explicit Models. Methods in Ecology and Evolution

# Here I assume that you already went to the quick-tutorial and the step-by-step-tutorial. If you didn't, go with those first.
  # any question: melisa.olave@uni-konstanz.de

# Remember to have the ms executable in the working directory.

# Lets run two different models using run.xl.pipeline in parallel. You will need foreach and doMC libraries. Make sure you have them installed.
library(foreach);
library(doMC);
library(ape);
library(gmp);

# set here your working directory and load functions:
setwd("/Users/Desktop/Olave_etal2017-MEE/examples");
source("functions.R");

# lets propose 2 models and set them into a data.frame:
  # model 0, no gene flow:
  model.number <- 0;
  time <- 0;
  species1 <- 0;
  species2 <- 0;
  M <- 0;
  model0 <- data.frame(model.number, time, species1, species2, M);

  # model 1, gene flow from B to A at time 0.5 and with M=1:
  model.number <- 1;
  time <- 0.5;
  species1 <- "B";
  species2 <- "A";
  M <- 1;
  model1 <- data.frame(model.number, time, species1, species2, M);

# now, put both models into a data.frame:
models <- rbind(model0, model1);
  # A suggestion: note that you could have an external txt file with a set of model description in each row (follow this "models" data frame as example), and load it using read.table function.

# Again... remember to have the ms executable in the working directory, and provide here the path to phylonet executable:
#PhylonetDir <- "myPath/phylonet/phylonet_v2_4.jar"

# run.xl.pipeline.parallel is identical to run.xl.pipeline (it will call it), except for models.table and cores
  # model.table = give here the table that we just created "models"
    # Note that col names in the model.table HAVE to be "model.number", "time", "species1", "species2" and "M".
  # cores = the number of cores you want to use to run in parallel. Each model will be ran in a different core, thus, because we have here 2 models, we will run in 2 cores:
source("/Users/melisaolave/Desktop/GitHub/Olave_etal2017-MEE/functions.R")


table <- run.xl.pipeline.parallel(wd=getwd(), model.table=models, cores=2, 
                          sptreeName="sp-tree.tre", genetreePattern="tre$",genetree.Path="genetrees", ImapName="Imap.txt",
                          brlength.correction=NULL, nloci=10, PhylonetDir=PhylonetDir,
                          saveHistogram=T, height=8.27, width=8.27);

# no progress can be reported using parallel cores! just be pacient!

# then check the table:
table

# Then you can performe the likelihood ratio test, get the likelihood of null and alternative models:
M0.likelihood <- table$likelihood[1];
M1.likelihood <- table$likelihood[2];

# no the likelihood ratio test, just 1 degree of freedom since there is only 1 extra parameter in model 1:
test <- Lratio.test(M0=M0.likelihood, M1=M1.likelihood, df=1);
test

# the model selected is the null mode M=0, since there is no great advance in likelihood calculations between models.


###### Final note:
# if want to include more than one M parameter, set the model.number to the same number, an the function will take the M in the same model.
  #For example, to simulate 2 migration parameters, one from A to B at time 0.5 and M=1, and the second from C to D at time 0.7 and M=2:
model.number <- c(2,2);
time <- c(0.5, 0.7);
species1 <- c("A", "C");
species2 <- c("B", "D");
M <- c(1,2);
model2 <- data.frame(model.number, time, species1, species2, M);

