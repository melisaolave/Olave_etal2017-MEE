# When using, please cite: Olave, M., Avila, L.J., Sites, J.W. Jr and Morando, M (2017). Detecting Hybridization by Likelihood Calculation of Gene Tree Extra Lineages given Explicit Models. Methodos in Ecology and Evolution
# This is an example of how to implement the method discribed by Olave et al. (2017; Methodos in Ecology and Evolution)
  # any question: melisa.olave@uni-konstanz.de

# REQUIREMENTS
  ##functions.R --> github.com/melisaolave/Olave_etal2017-MEE
  ##ape library installed in R
  ##ms (Hudson 2002) located at the same folder than the data
  ##Phylonet v2.3 (Than and Nakhleh 2009) # will not work with a different version!

# Previous files you should have ready
  ## a species tree in newick format and with branch length in coalescent units. See example sp-tree.tre
  ## an Imap file containing the individual-species asociations. See Imap.txt as example.

# The sequence of steps from (i) to (iii) is shown in (Fig. 2; Olave et al. in press) and corresponds to:
# (i) Proposing a model including a species tree topology and branch length in coalescent units (CU = t/Ne; where t is generations and Ne is the effective population size). Tip: Species tree programs like ASTRAL (Mirarab et al. 2014) and MP-EST (Liu et al. 2010) will now calculate the species tree in CU.
# (ii) Inference of expected extra lineages given the model. A set of H gene trees is simulated following the model described in (i). Extra lineages are counted between each simulated gene tree and the species tree S proposed in (i), providing a distribution of the expected number of extra lineages under the model in (i). The method currently uses the software ms (Hudson 2002). Then the package Phylonet (Than and Nakhleh 2009) is used to count extra lineages, and provide the distribution of expected discordance.
# (iii) Likelihood of the empirical data. The empirical gene trees are compared to the distribution generated in (ii) and the likelihood calculations are obtained.
    #Steps i - iii should be repeated with as many models as is desired.
# (iv) Model selection. A likelihood ratio test is then used to select the best model among the proposed ones.

# make sure you have installed the ape library 
library(ape);

# set your working directory
setwd("/User/Desktop/Olave_etal2017-MEE/examples");
wd <- getwd(); # getting rooth directory

# load functions
source("functions.R");

# Follwing steps i - iv are described in Olave et al. in press (see figure 2) and README.txt

# Step (i) 
## get --> species tree topology in newick format with branch length in coalescent units.
tree <- read.tree("sp-tree.tre");

# get the species names
species <- tree$tip.label;
species <- species[order(species)]; # it is important to sort them for getting the right command below

# getting total number of species
nspecies <- length(species);

# read the sp tree as text, to later use it at input for get.topo2ms function
tree.newick <- scan("sp-tree.tre", what="character", sep="\n");

## change newick format to ms program format
topology <- get.topo2ms(tree.newick, species, brlength.correction=NULL);

# NOTE: always chech if the topology is in the right format. For this case, it should look like this: "-ej 1 2 1 -ej 2 3 1 -ej 3 4 1"
topology

# provide the M parameter: this case M = 0, so we can get the likelihood of the null model (i.e. no gene flow)
ms.gf.code <- NULL;

# Step (ii). Generating expectations given the model in (i)
# getting Imap
Imap.txt <- read.table("Imap.txt", sep="\t", header=T);

# get Imap in required format (with individuals named as numbers)
Imap <- get.Imap(indVec, nspecies);

# getting number of individuals per species into indVec vector - in this case all species has 3 individuals, but you could have different numbers, so this loops will give you the sampling number for each species
indVec <- NULL
for(i in 1:length(species)){
  subImap <- Imap.txt[Imap.txt$species == species[i],];
  indVec[i] <- nrow(subImap);
}

# getting taxaVec
totalTaxa <- sum(indVec);
taxaVec <- paste(indVec, collapse=" ");


# define H
nloci <- 10000; # number of loci or gene trees that want to simulate (H). Recommended = 10,000

# name your ms output file (will contain all simulated gene trees for expectations). It is usually recommended to add as much details as possible of the simulated parameters: number of individuals per species and M paramenter (in this M = 0);
Exp.outputName <- paste("expectedGeneTrees-", paste(species, indVec,  collapse="", sep=""), # this part will collapse the species names with the number of individuals per each of them
                        "-M0.txt", sep=""); # M = 0

# executing ms. After this, you should get a file named as it was defined in Exp.outputName above. Check it to see if it has right gene trees
ms(totalTaxa=totalTaxa, nloci=nloci, species=species, taxaVec=taxaVec, ms.topology=topology,
   gf.code=ms.gf.code, outputName=Exp.outputName);

# you could have the program PhyloNet at a different path. If so, name it here by replacing wd below
PhylonetDir <- file.path(wd, "phylonet_v2_3.jar"); # only phylonet_v2_3.jar version accepted!

# Time to get the number of extra lineages observed!
distVec <- sptree.vs.genetrees(wd=wd, PhylonetDir=PhylonetDir, tree.newick=tree.newick, 
                               genetreeName=Exp.outputName, species=species, Imap=Imap, sp.treeName="sp-tree.tre");

# if you want, you can save the progress by saving what is recorded in distVec (then you could read it laterif want to repeat analyses and save all this time so far);
distVecName <- paste("distVec-",Exp.outputName, sep="");
write(distVec, distVecName);

# NOTE: If you ever want to load the distVec information from the file do this:
distVecTable <- read.table(distVecName);
for(i in 1:ncol(distVecTable)){
  distVec <- unlist(c(distVec, distVecTable[,i])); 
}

# you can check the distribution of expected extra lineages by plotting distVec
hist(distVec,  xlab = 'expected extra lineages', ylab = "Frequency", main=NULL);

# Step iii - getting likelihood of the real gene trees

# get your real gene trees - idealy they will be at a folder "genetrees" at the same directory
setwd(file.path(wd, "genetrees"));
my.genetrees <- "tre" # provide a pattern that represent all names, for example, all the example gene trees has "tre"
genetreeFiles <- list.files(pattern = my.genetrees); # this will get the gene tree names

# make sure that you got the right gene trees (and that the species tree is not in the list!)
genetreeFiles

# get likelihood
  # In case you have tried this several times, make sure that you are providing the right distVec here!!!
setwd(wd);
likelihood <- getLikelihood(sptree = new.sptree, genetreeFile = genetreeFiles, 
                            distVec=distVec, 
                            sp.treeName = sp.treeName,
                            PhylonetDir = PhylonetDir, wd = wd,
                            Imap = Imap.txt, species = species, taxaVec = taxaVec);
# you will see the histogram of expected extra lineages given the model in step i and comparision with empirical observation (your data = red line)

# the log likelihood of this model is
likelihood

# get the information so far into a table
result <- data.frame(distVecName, "M=0", likelihood);
colnames(result) <- c("distVec.Source", "M", "likelihood");

# NOTE: the steps so far need to be repeated with a different model, including some magnitutes of M, for example try:
# When having M different to 0, you can use gf.code2ms to automathically convert the format
  # time = time when gene flow occured back in time (then 0 means current gene flow);
  # species1 = from
  # species2 = recieving species
  # M= migration parameter magnitute
ms.gf.code <- gf.code2ms(time=0, species1="A", species2="B", M=2, all.species=species);

# if you want to add one more M, for example from B to A, with M=1
ms.gf.code[2] <- gf.code2ms(time=0, species1="B", species2="A", M=1, all.species=species);

# remember to change the output name with the new M information:
Exp.outputName <- paste("expectedGeneTrees-", paste(species, indVec,  collapse="", sep=""), # this part will collapse the species names with the number of individuals per each of them
                        "-M_AB2-BA1",".txt", sep="");

# then run ms to simulated the expected gene trees under this new model
ms(totalTaxa=totalTaxa, nloci=nloci, species=species, taxaVec=taxaVec, ms.topology=topology,
   gf.code=ms.gf.code, outputName=Exp.outputName);

# get extra lineages
distVec <- sptree.vs.genetrees(wd=wd, PhylonetDir=PhylonetDir, tree.newick=tree.newick, 
                               genetreeName=Exp.outputName, # always make sure that you are getting the right ms simulations
                               species=species, Imap=Imap, sp.treeName="sp-tree.tre");
distVecName <- paste("distVec-",Exp.outputName, sep=""); # save distVec in case that you want to load it later to repeat calculations
write(distVec, distVecName);

# get likelihood for this model
likelihood <- getLikelihood(sptree = new.sptree, genetreeFile = genetreeFiles, 
                            distVec=distVec, 
                            sp.treeName = sp.treeName,
                            PhylonetDir = PhylonetDir, wd = wd,
                            Imap = Imap.txt, species = species, taxaVec = taxaVec);

# add the new information into the previous table
new.results <- c(distVecName, "M=0", likelihood);
result <- rbind(result, new.results);
result

# Lets say that these are all the models that we are interested to test, thus you can perfome a likelihood ratio test between this two models
  # the df are the degrees of freedom, and this depends in how many M parameters you have. Here we included 2 M parameters (gene flow from A to B and from B to A)
test <- Lratio.test(M0=result$likelihood[1], M1=result$likelihood[2], df=2);
test

# the result is a list with the G statistic calculation, the p-value of the hypothesis M0 = M1, and the model selected (if M=0 or M>0);















