# When using, please cite: Olave, M., Avila, L.J., Sites, J.W. Jr and Morando, M (2017). Detecting Hybridization by Likelihood Calculation of Gene Tree Extra Lineages given Explicit Models. Methods in Ecology and Evolution

# This is a quick tutorial example of how to implement the method discribed by Olave et al. (2017; Methods in Ecology and Evolution)
  # Here, I automated all the pipeline in a single function. However, note that there are many steps that can be done per separate, and it is not a bad idea to check the extended tutorial step-by-step-tutorial
  # any question: melisa.olave@uni-konstanz.de

# REQUIREMENTS
  ##functions.R --> github.com/melisaolave/Olave_etal2017-MEE
  ##ape and gmp libraries installed in R
  ##ms (Hudson 2002) compiled and located at the same folder than the data
  ##Phylonet v2.4 (Than and Nakhleh 2009) # it might not work with a different version!

# Files you should have ready
  ## a species tree in newick format and with branch length in coalescent units. See example sp-tree.tre
  ## an Imap file containing the individual-species associations. See Imap.txt as example.
  ## a set of gene trees. See examples in genetree1.tre - genetree5.tre 

# This quick tutorial allows to automatise all the sequence of steps from (i) to (iii) is shown in (Fig. 2; Olave et al. 2017 - MEE) and corresponds to:
  # (i) Proposing a model including a species tree topology and branch length in coalescent units (CU = t/Ne; where t is generations and Ne is the effective population size). Tip: Species tree programs like ASTRAL (Mirarab et al. 2014) and MP-EST (Liu et al. 2010) will now calculate the species tree in CU.
  # (ii) Inference of expected extra lineages given the model. A set of H gene trees is simulated following the model described in (i). Extra lineages are counted between each simulated gene tree and the species tree S proposed in (i), providing a distribution of the expected number of extra lineages under the model in (i). The method currently uses the software ms (Hudson 2002). Then the package Phylonet (Than and Nakhleh 2009) is used to count extra lineages, and provide the distribution of expected discordance.
  # (iii) Likelihood of the empirical data. The empirical gene trees are compared to the distribution generated in (ii) and the likelihood calculations are obtained.
  #Steps i - iii should be repeated with as many models as is desired.
  # (iv) Model selection. A likelihood ratio test is then used to select the best model among the proposed ones.


#### overview to the run.xl.pipeline function:
# wd=getwd()      ## working directory path
# sptreeName="sp-tree.tre"      ## species tree name
# genetreePattern="tre$"      ## a pattern to match gene tree files (with reg expressions)
# genetree.Path=getwd()      ## a path to gene trees files
# ImapName="Imap.txt"      ## the Imap name, indicating individual species associations. Columns should be named "traits" and "species", see Imap.txt example.
# brlength.correction=NULL      ## in case the branch length in the species tree are not in coalescent units, or if gene trees come from a mitochondrial or chloroplast locus, you can correct the branch length here
# nloci=100      ## number of loci to simulate with ms and produce the distribution of expected extra lineages. Default (and recommended)= 10000;
# PhylonetDir="phylonet_v2_3.jar"      ## a path where the Phylonet directory and executable are
# PhyloImapName="PhyloImap.txt"      ## the name of the PhyloImap output that will be saved during the pipeline
# geneTreeTempName="genetree.temp"     ## the name for a temporal file that will be created during the pipeline
# saveXLVec=T     ## True if want to save the extra lineage vector estimated to infer the distribution. Recommended.
# loadXLvec=F     ## if True, will load a previous saved extra lineage vector for the given model
# only.XLvec=F    ## if True, then only the extra lineage vector is estimated, and no comparision with true gene trees is performed.
# time=0     ## vector of time of hybridization event. 0 for no hybridization.
# species1=0     ## vector of species names for hybridization event (giver). 0 for no hybridization.
# species2=0     ## vector of species names for hybridization event (recipient). 0 for no hybridization.
# M=0     ## vector of migration parameter magnitude. 0 for no hybridization.
# saveHistogram=T     ## True if want to save a plot of the histogram of extra lineage given the model.
# height=8.27      ## height of the histogram
# width=height     ##width of the histogram


####################### Below a quick tutorial  #########################
# make sure you have installed the ape and gmp libraries 
library(ape);
library(gmp);

# Download the functions.R and load them:
setwd("/Users/Desktop/Olave_etal2017-MEE/");
source("functions.R");

# provide Phylonet executable:
PhylonetDir <- "mypath/phylonet/phylonet_v2_3.jar";

# once you have the species tree, gene trees, Imap and ms executable in the working directory, lets get the likelihood for a null model (no gene flow):
  # for a null model here we are setting time=0, species1=0, species2=0, M=0
table <- run.xl.pipeline(wd=getwd(), sptreeName="sp-tree.tre", genetreePattern="tre$",genetree.Path="genetrees", ImapName="Imap.txt",
                      brlength.correction=NULL, nloci=50, PhylonetDir=PhylonetDir,
                      PhyloImapName="PhyloImap.txt", geneTreeTempName="genetree.temp", saveXLVec=T, loadXLvec=F, only.XLvec=F,
                      time=0, species1=0, species2=0, M=0,
                      saveHistogram=T, height=8.27, width=8.27);

# you should see two new files in your working directory folder:
  # 1. the file XLvec-XLexpected-A3B3C3D3-M0.txt has the extra lineage vector required to reconciliate for the total of 100 loci simulated to infer the distribution (step ii, described above).
    #Note that the name of the file includes de number of individuals per species detected A3,B3,C3,D3 (in this case, all species have 3 individuals),
    # and the model description, in this case no gene flow, so M0. If there is gene flow you will see a description of the model, including time, species and M parameter.
  # 2. the Histogram-XLexpected-A3B3C3D3-M0.txt.pdf is an histogram constructed with the extra lineage vector

# also, the table has now the likelihood of this null model, check it:
table

# now you can a run a new model, including gene flow, from A to B species in time 0.5 and M=1, and from C to D species in time 0.7 and M=2:
  # lets difine:
time.vec <- c(0.5,0.7);
species1.vec <- c("A", "C");
species2.vec <- c("B", "D");
M.vec <- c(1,2);
 #  then run:
new.table <- run.xl.pipeline(wd=getwd(), sptreeName="sp-tree.tre", genetreePattern="tre$", genetree.Path="genetrees", ImapName="Imap.txt",
                         brlength.correction=NULL, nloci=100, PhylonetDir=PhylonetDir,
                         PhyloImapName="PhyloImap.txt", geneTreeTempName="genetree.temp", saveXLVec=T, loadXLvec=F,
                         time=time.vec, species1=species1.vec, species2=species2.vec, M=M.vec,
                         saveHistogram=T, height=8.27, width=height);

# now should get two new files (as before, but corresponding to the new model), as well as the new table:
new.table

# it is posible to re run the run.xl.pipeline functions with as many different models and parameter magnitudes as desire.
  # then you can conviene all results:
final.table <- rbind(table, new.table);

# Then you can performe the likelihood ratio test, get the likelihood of null and alternative models:
M0.likelihood <- final.table$likelihood[1];
M1.likelihood <- final.table$likelihood[2];

# in this case the alternative model has 2 M parameters (gene (1) from A to B and (2) from C to D), thus the difference of parameter numbers between M0 and M1 is = 2. These are the degrees of freedom:
test <- Lratio.test(M0=M0.likelihood, M1=M1.likelihood, df=2);
test

# the model selected is M=0, no gene flow detected for this dataset.

