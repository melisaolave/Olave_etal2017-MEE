# get.Imap function
get.Imap <- function(fixed.nsamples = 5, nspecies = 8){
  if(length(fixed.nsamples) == 1){
    traits <- 1:(nspecies*fixed.nsamples);
    Imap.species <- character();
    species <- LETTERS[1:nspecies]
    for(i in species){
      for (j in 1:fixed.nsamples){
        Imap.species <- c(Imap.species, i);
      }
    }
  }
  if(length(fixed.nsamples) > 1){
    traits <- integer();
    Imap.species <- character();
    species <- LETTERS[1:nspecies]
    count <- 0;
    count2 <- 1;
    for(k in fixed.nsamples){
      traits1 <- (count+1):(k+count);
      traits <- c(traits, traits1);
      count <- count + k
      for (i in 1:k){
        Imap.species <- c(Imap.species, species[count2]);
      }
      count2 <- count2 +1
    }
  }
  Imap <- data.frame(traits, Imap.species);
  colnames(Imap) <- c("traits", "species");
  return(Imap)
}

# get topological info from tree and codify to ms

get.topo2ms <- function(tree.newick, species, brlength.correction = NULL){
  count <- 1;
  new.tree.newick <- tree.newick;
  for(i in species){ # to replace species names by numbers, then ms can take them
    new.tree.newick <- gsub(pattern=i, replacement=count, x=new.tree.newick);
    count <- count+1;
  }
  topology <- character();
  for(j in 1:length(tree.newick)){ # in case you provided more than one sp tree topology
    if(grepl(pattern=";", x=new.tree.newick[j]) == F){ # if your sp tree doesn't have a ";" at the end, then add it
      new.tree.newick[j] <- gsub(pattern="(.*)$", replacement="\\1;", x=new.tree.newick[j]);
    }
    taxa1.Vec <- character();
    taxa2.Vec <- character();
    blength1.Vec <- numeric();
    blength2.Vec <- numeric();
    i <- T;
    while(i){
      if(length(unlist(strsplit(gsub(pattern="[^,]", replacement="", new.tree.newick), ","))) == 1){
        clade <- gsub(pattern="^.*(\\([[:alnum:]]+:\\d\\.*\\d*,[[:alnum:]]+:\\d\\.*\\d*\\):\\d\\.*\\d*).*", 
                      replacement="\\1;", x=new.tree.newick[j]);
      }else{
        clade <- gsub(pattern="^.*(\\([[:alnum:]]+:\\d\\.*\\d*,[[:alnum:]]+:\\d\\.*\\d*\\):\\d\\.*\\d*).*", 
                      replacement="\\1", x=new.tree.newick[j]);
      }
      if(length(unlist(strsplit(gsub(pattern="[^,]", replacement="", clade), ","))) > 1){ # if there is more than one "," in clade, means a polytomy
        #clade <- gsub(pattern="")
        cladeBackup <- clade;
        while(length(unlist(strsplit(gsub(pattern="[^,]", replacement="", clade), ","))) > 1){
          taxa1 <- gsub(pattern="^\\(([[:alnum:]]+):.*", replacement="\\1", x=clade);
          taxa2 <- gsub(pattern="^\\(.*,([[:alnum:]]+):.*", replacement="\\1", x=clade);
          blength1 <- as.numeric(gsub(pattern=paste("^\\(.*,[[:alnum:]]+:.*", taxa2, ":(\\d+\\.*\\d*).*", sep=""), replacement="\\1", x=clade));
          
          taxa1.Vec <- c(taxa1.Vec, taxa1);
          taxa2.Vec <- c(taxa2.Vec, taxa2);
          blength1.Vec <- c(blength1.Vec, blength1);
          clade <- gsub(pattern=paste("^(\\(", taxa1, ":.*),", taxa2, ":.*(\\).*)", sep=""),#"^.*(\\([[:alnum:]]+:\\d\\.*\\d*,[[:alnum:]]+:\\d\\.*\\d*\\):\\d\\.*\\d*).*", 
                        replacement="\\1\\2", x=clade);     
        }
        new.clade <- clade;
        cladeBackup <- gsub(pattern="\\(", replacement="\\\\(", cladeBackup);
        cladeBackup <- gsub(pattern="\\)", replacement="\\\\)", cladeBackup);
        new.tree.newick[j] <- gsub(pattern=paste("(^.*)(", cladeBackup, ")(.*)", sep=""),            
                                   replacement=paste("\\1", new.clade, "\\3", sep=""), x=new.tree.newick[j]);
      }else if(grepl(pattern=";", x=clade) == F){
        taxa1 <- gsub(pattern="^\\(([[:alnum:]]+):.*", replacement="\\1", x=clade);
        taxa2 <- gsub(pattern="^\\(.*,([[:alnum:]]+):.*", replacement="\\1", x=clade);
        blength1 <- as.numeric(gsub(pattern="^\\([[:alnum:]]+:(.*),.*", replacement="\\1", x=clade));
        taxa1.Vec <- c(taxa1.Vec, taxa1);
        taxa2.Vec <- c(taxa2.Vec, taxa2);
        blength1.Vec <- c(blength1.Vec, blength1);
        blength2 <- as.numeric(gsub(pattern="^\\(.*,[[:alnum:]]+:.*:(\\d+\\.*\\d*)", replacement="\\1", x=clade));
        blength2.Vec <- c(blength2.Vec, blength2);
        new.blength <- blength1 + blength2;
        new.clade <- paste(taxa1, ":", new.blength, sep="");
        new.tree.newick[j] <- gsub(pattern="(^.*)(\\([[:alnum:]]+:\\d\\.*\\d*,[[:alnum:]]+:\\d\\.*\\d*\\):\\d\\.*\\d*)(.*)",            
                                   replacement=paste("\\1", new.clade, "\\3", sep=""), x=new.tree.newick[j]);
      }else{
        taxa1 <- gsub(pattern="^\\(([[:alnum:]]+):.*", replacement="\\1", x=clade);
        taxa2 <- gsub(pattern="^\\(.*,([[:alnum:]]+):.*", replacement="\\1", x=clade);
        blength1 <- as.numeric(gsub(pattern="^\\([[:alnum:]]+:(.*),.*", replacement="\\1", x=clade));
        taxa1.Vec <- c(taxa1.Vec, taxa1);
        taxa2.Vec <- c(taxa2.Vec, taxa2);
        blength1.Vec <- c(blength1.Vec, blength1);
        new.blength <- blength1;
        new.clade <- paste(taxa1, ":", new.blength, sep="");
        new.tree.newick[j] <- gsub(pattern="(^.*)(\\([[:alnum:]]+:\\d\\.*\\d*,[[:alnum:]]+:\\d\\.*\\d*\\):\\d\\.*\\d*)(.*)",            
                                   replacement=paste("\\1", new.clade, "\\3", sep=""), x=new.tree.newick[j]);
        i <- F;
      }
    }
    if(length(brlength.correction) != 0){ # this is to correct branch length unites if sp tree was estimated by a phylo inference program
      blength1.Vec <- blength1.Vec/brlength.correction;
    }
    blength1.Vec <- round(blength1.Vec, 4);
    topology[j] <- paste("-ej", blength1.Vec, taxa2.Vec, taxa1.Vec, collapse = " ");
  }
  return(topology)
}


# gene flow and pop size changes to ms

gf.code2ms <- function(time=0, species1, species2, M=0, all.species){
  # -em t i j x
    # t = time; i = species i; j = species j; x = migration parameter (4Nm)
  gf.code <- paste(time, species1, species2, M);
  count <- 1;
  gf.code <- gsub(pattern="(\\w+)$", replacement="\\1 ", x=gf.code); # to add an space at the end of gf.code provided
  for(i in all.species){ # to replace species names by numbers, then ms can take them
    gf.code <- gsub(pattern=paste("([[:space:]])",i, "([[:space:]])", sep=""), # the space is just to be sure it's matching taxa names, even when taxa names are numbers.
                    replacement=paste("\\1",count,"\\2", sep=""), 
                    x=gf.code);
    count <- count+1;
  }
  gf.code <- gsub(pattern="(\\w+)[[:space:]]$", replacement="\\1", x=gf.code); # to delete space at the end of code
  gf.code <- paste("-em", gf.code);
  return(gf.code);
}

# execute ms
ms <- function(totalTaxa, nloci, species, taxaVec, ms.topology, gf.code = NULL, popSize.code = NULL,  ms.more.arg = NULL, outputName = "genetrees.txt", ms.exec="./ms", rm.header=T){
  if(length(gf.code) >1){
    gf.code <- paste(gf.code, collapse=" ");
  }
  if("ms" %in% list.files() == F){
    stop("ms executable not found in working directory.
         1. Download ms program: http://home.uchicago.edu/rhudson1/source/mksamples.html
         2. Compile it.
         3. Paste ms executable in working directory:", wd, "\n...and run again!\n");
  }
  command <- paste(ms.exec, 
                   totalTaxa, 
                   nloci, 
                   "-T -I", 
                   length(species), 
                   taxaVec, 
                   ms.topology,
                   gf.code,
                   popSize.code,
                   ms.more.arg,
                   ">", 
                   outputName);
  cat("ms command:", command, "\n");
  system(command);
  if(rm.header){
    ms.out <- scan(outputName, what="character", sep="\n");
    ms.out <- ms.out[4:length(ms.out)];
    ms.out <- gsub(pattern="//", replacement="", ms.out);
    write(ms.out, outputName);
  }
}

### write PhyloImap
writePhyloImap <- function(Imap, species, msTaxaNames = F, taxaVec = NULL, PhyloImapName="PhyloImap.txt"){
  PhyloImap <- character();
    # writting Imap required by Phylonet
  if(msTaxaNames){
    if(length(taxaVec) == 0){
      stop("To write PhyloImap using ms taxa names, you need to provide a taxaVec\n");
    }
    if(is.integer(taxaVec) == FALSE){
      taxaVec <- as.integer(unlist(strsplit(taxaVec, split=" ")));
    }
    count <- 1;
    start <- 1;
    for(x in taxaVec){
      PhyloImap[count] <- paste(species[count], ":", paste(start:(start-1+x), collapse=","), ";", sep="");
      start <- start + x;
      count <- count +1;
    }
  }else{
    for(z in 1:length(species)){
      taxa <- Imap[Imap$species == species[z],];
      taxa <- taxa$traits;
      taxa <- paste(taxa, collapse=",");
      PhyloImap[z] <- paste(species[z], ":", taxa, ";", sep="");
    }
  }
  write(PhyloImap[1:length(species)], PhyloImapName);
}

#### sp tree vs genetrees
sptree.vs.genetrees <- function(wd= getwd(), PhylonetDir = getwd(), PhyloImapName="PhyloImap.txt",
                               tree.newick, genetreeNames, species, Imap, msTaxaNames = F, taxaVec = NULL, sp.treeName, sptree.dir = getwd(), geneTreeTempName="genetree.temp"){
  writePhyloImap(Imap=Imap, species=species, msTaxaNames = msTaxaNames, taxaVec = taxaVec, PhyloImapName = PhyloImapName);
  distVec <- integer();
  cat("Reading expected gene trees from file:",genetreeNames,"\n");
  genetrees <- read.tree(genetreeNames);
  if(class(genetrees) == "multiPhylo"){
    cat("Estimating expected distribution of extra lineages, be patient...\n");
    count <- 1;
    for(l in 1:length(genetrees)){
      write.tree(genetrees[[l]], geneTreeTempName);
      command <- paste("java -jar ", PhylonetDir, " deep_coal_count ", sptree.dir, "/", sp.treeName," ", wd, "/", geneTreeTempName, " -a ", wd, "/", PhyloImapName, sep="");
      deepcoal <- system(command, intern = T);
      deepcoal <- as.integer(gsub(pattern=".*lineages: (\\d+)", replacement="\\1", deepcoal[2]));
      distVec <- c(distVec, deepcoal);
      cat("progress: ",l , " / ", length(genetrees), "\n");
      count <- count +1;
      file.remove(geneTreeTempName);
    }
  }else{
    command <- paste("java -jar ", PhylonetDir, " deep_coal_count ", wd, "/", sp.treeName," ", wd, "/", genetreeNames, " -a ", wd, "/", PhyloImapName, sep="");
    deepcoal <- system(command, intern = T);
    deepcoal <- as.integer(gsub(pattern=".*lineages: (\\d+)", replacement="\\1", deepcoal[2]));
    distVec <- c(distVec, deepcoal);
  }
  file.remove(PhyloImapName);
  return(distVec)
}

get.myGenetreeXL <-function(genetreeFiles, genetree.Path=getwd(), sptree.dir=getwd(), sp.treeName="sp-tree.tre",species, Imap="Imap.txt", 
                            msTaxaNames=F, taxaVec=NULL, PhyloImapName="PhyloImap.txt", PhylonetDir=getwd()){
  myDist <- NULL;
  for(i in 1:length(genetreeFiles)){
    mygenetreeFile <- file.path(genetree.Path, genetreeFiles[i]);
    writePhyloImap(Imap, species, msTaxaNames = F, taxaVec = F, PhyloImapName=PhyloImapName);
    # comparing with real data
    command <- paste("java -jar ", PhylonetDir, " deep_coal_count ", sptree.dir, "/", sp.treeName," ", mygenetreeFile, " -a ", PhyloImapName, sep="");
    mydeepcoal <- system(command, intern = T);
    myDist <- c(myDist, as.integer(gsub(pattern=".*lineages: (\\d+)", replacement="\\1", mydeepcoal[2])));
  }
  return(myDist);
}
    
# comparing real gene trees with model
getLikelihood <- function(genetreeFiles, distVec, sp.treeName = NULL,
                          PhylonetDir = getwd(), wd = getwd(), genetree.Path=getwd(),
                          Imap, species, msTaxaNames=F, taxaVec=NULL, sptree.dir=getwd(), PhyloImapName="PhyloImap.txt"){
  library(gmp);
  pVec <- numeric();
  nVec <- integer();
  myDistVec <- integer();
  H <- length(distVec);
  N <- length(genetreeFiles);
  for(i in 1:N){
    if(class(genetreeFiles) == "character"){
      myDist <- get.myGenetreeXL(genetreeFiles=genetreeFiles[i], genetree.Path=genetree.Path, sptree.dir=sptree.dir, sp.treeName=sp.treeName,species=species, Imap=Imap, 
                                  msTaxaNames=msTaxaNames, taxaVec=taxaVec, PhyloImapName=PhyloImapName, PhylonetDir=PhylonetDir);
    }else{
      myDist <- genetreeFiles[i];
    }
    cat("Extra lineages observed:", myDist, "\n");
    if(all(myDist != myDistVec)){
      h <- sum(distVec == myDist);
      if(h == 0){ # this is to prevent p = 0, then we choose a very small number
        pVec <- c(pVec, 0.0000000001);
      }else{
        p <-  h/H;
        hist(distVec,  xlab = 'expected extra lineages', ylab = "Frequency", main=NULL);
        abline(v = myDist, col = "red");
        pVec <- c(pVec, p); 
      }
    }
    myDistVec[i] <- myDist;
  }
  my.freq <- unique(myDistVec);
  count <- 1;
  for(j in my.freq){
    nVec[count] <- sum(myDistVec == j);
    count <- count+1;
  }
  file.remove(PhyloImapName);
  # N!/ni! * pi^ni
    # ln N! - ln ni! + ni * ln pi
  likelihood <- log(factorialZ(N)) - sum(log(factorialZ(nVec))) + sum(nVec*log(pVec));
  return(likelihood);
}

## Likelihood ratio test between two given models
Lratio.test <- function(M0, M1, df){
  M0 <- as.numeric(M0);
  M1 <- as.numeric(M1);
  G <- 2*(M1 - M0);
  ratio.test <- pchisq(G, df, lower.tail = F);
  if(ratio.test >= 0.05){
    model.selected <- "M=0"; 
    bestLikelihood <- M0;
  }else{ 
    model.selected <- "M>0";
    bestLikelihood <- M1;
  }
  result <- list(G, ratio.test, bestLikelihood, model.selected);
  names(result) <- c("G", "p-value", "likelihood","model.selected");
  return(result);
}

# get individual number from a fasta matrix and codify ms code
# this function can be useful if you want to increase accurancy because of missing data, but this is not implemented in the tutorial provided
getTaxa2ms <- function(fileName, species, Imap){
  # read fasta matrix
  matrix <- scan(fileName, what="character", sep="\n");
  seq.starts <- grep(pattern=">", x=matrix);
  seq <- character();
  taxa <- character();
  taxaVec <- integer();
  for(j in 1:length(seq.starts)){ # to get taxa names vector
    if(j != max(length(seq.starts))){
      sequence <- matrix[(seq.starts[j]+1):(seq.starts[j+1]-1)];
    }else{
      sequence <- matrix[(seq.starts[j]+1):length(matrix)];
    }
    taxa[j] <- matrix[seq.starts[j]];
  }
  taxa <- gsub(pattern=">", replacement="", x=taxa); # deleting ">"
  
  for (k in species){ # to recognize number of individuals to be simulated per species
    subset <- Imap[Imap$species == k,];
    ind <- as.character(subset$traits);
    matched.taxa <- sum(ind %in% taxa);
    taxaVec <- c(taxaVec, matched.taxa); # vector of number of individual per species from seq matrices
  }
  return <- taxaVec
}

### get.pvalue  - this function will return only the probability vector (could be needed when having missing data)
get.p <- function(sptree, genetreeFile, distVec, sp.treeName = NULL,
                  PhylonetDir = getwd(), sp.tree.wd = getwd(), genetree.wd = getwd(),
                  Imap, species, taxaVec, PhyloImapName="PhyloImap.txt"){
  pVec <- numeric();
  nVec <- integer();
  myDistVec <- integer();
  H <- length(distVec);
  N <- length(genetreeFile);
  for(i in 1:N){
    mygenetreeFile <- genetreeFile[i];
    writePhyloImap(Imap, species, msTaxaNames = F, taxaVec = NULL, PhyloImapName=PhyloImapName);
    # comparing with real data
    command <- paste("java -jar ", PhylonetDir, " deep_coal_count ", sp.tree.wd, "/", sp.treeName," ", genetree.wd, "/", mygenetreeFile, " -a ", PhyloImapName, sep="");
    mydeepcoal <- system(command, intern = T);
    file.remove(PhyloImapName);
    myDist <- as.integer(gsub(pattern=".*lineages: (\\d+)", replacement="\\1", mydeepcoal[2]));
    cat("Extra lineages observed:", myDist, "\n");
    if(all(myDist != myDistVec)){
      h <- sum(distVec == myDist);
      if(h == 0){ # this is to prevent p = 0, then we choose a very small number
        pVec <- c(pVec, 0.0000000001);
      }else{
        p <-  h/H;
        hist(distVec,  xlab = 'expected extra lineages', ylab = "Frequency", main=NULL);
        abline(v = myDist, col = "red");
        pVec <- c(pVec, p); 
      }
    }
    myDistVec[i] <- myDist;
  }
  output <- list(pVec,myDistVec);
  names(output) <- c("pVec", "myDistVec");
  return(output);
}

# get likelihood when providing a pVec, nVec and N
get.L <- function(nVec, pVec, N){
  # N!/ni! * pi^ni
  # ln N! - ln ni! + ni * ln pi
  likelihood <- log(factorialZ(N)) - sum(log(factorialZ(nVec))) + sum(nVec*log(pVec));
  return(likelihood);
}

getTaxaVec <- function(mygenetrees, Imap){
  indList <- list()
  Imap <- as.data.frame(Imap);
  species <- as.character(unique(Imap$species));
  for(i in 1:length(mygenetrees)){
    genetree <- read.tree(mygenetrees[i]);
    tips <- sort(genetree$tip.label);
    new.Imap <- Imap[as.character(Imap$traits) %in% tips,];
    indVec <- NULL;
    for(j in 1:length(species)){
      indVec[j]<- sum(new.Imap$species == species[j]);
    }
    indList[[i]] <- indVec;
  }
  return(indList);
}

get.n <- function(myDistVec){
  my.freq <- unique(myDistVec);
  count <- 1;
  nVec <- integer();
  for(j in my.freq){
    nVec[count] <- sum(myDistVec == j);
    count <- count+1;
  }
  return(nVec);
}

run.xl.pipeline <- function(wd=getwd(), sptreeName, genetreePattern="tre$",genetree.Path=getwd(), ImapName="Imap.txt",
                            ms.exec="./ms", brlength.correction=NULL, nloci=10000, PhylonetDir=NULL,
                            PhyloImapName="PhyloImap.txt", geneTreeTempName="genetree.temp", saveXLVec=T, loadXLvec=F, only.XLvec=F,
                            time=0, species1=0, species2=0, M=0,
                            saveHistogram=T, height=8.27, width=8.27){
  library(ape);
  library(gmp);
  if(length(time) != length(species1) | length(species2) != length(M) | length(time) != length(M)){
    stop("The vectors time, species1, species2 and M should have the same length\n");
  }
  setwd(wd);
  tree <- read.tree(sptreeName);
  species <- tree$tip.label;
  species <- species[order(species)]; 
  nspecies <- length(species);
  tree.newick <- scan(sptreeName, what="character", sep="\n");
  topology <- get.topo2ms(tree.newick, species, brlength.correction=brlength.correction);
  cat("topology:", topology, "\n")
  Imap.txt <- read.table(ImapName, sep="\t", header=T);
  indVec <- NULL
  for(i in 1:length(species)){
    subImap <- Imap.txt[Imap.txt$species == species[i],];
    indVec[i] <- nrow(subImap);
  }
  Imap <- get.Imap(indVec, nspecies);
  totalTaxa <- sum(indVec);
  taxaVec <- paste(indVec, collapse=" ");
  ms.gf.code <- NULL;
  if(any(M != 0)){
    for(i in 1:length(time)){
      if(M[i] != 0){
        ms.gf.code <- c(ms.gf.code, gf.code2ms(time=time[i], species1=species1[i], species2=species2[i], M=M[i], all.species=species));
      }
    }
  }  
  if(length(ms.gf.code) == 0){
    gfDesc <- "M0";
  }else{
    gfDesc <- paste(species1, species2,"_t=", time, "_M=", M, sep="", collapse="+");
  }
  Exp.outputName <- paste("XLexpected-", paste(species, indVec,  collapse="", sep=""), # this part will collapse the species names with the number of individuals per each of them
                          "-", gfDesc,".txt", sep=""); # M = 0
  distVecName <- paste("XLvec-",Exp.outputName, sep="");
  if(loadXLvec == F){
    ms(totalTaxa=totalTaxa, nloci=nloci, species=species, taxaVec=taxaVec, ms.topology=topology,
       gf.code=ms.gf.code, outputName=Exp.outputName, ms.exec=ms.exec);
    distVec <- sptree.vs.genetrees(wd=wd, PhylonetDir=PhylonetDir, PhyloImapName=PhyloImapName, tree.newick=tree.newick, 
                                   genetreeNames=Exp.outputName, species=species, Imap=Imap, sp.treeName=sptreeName, geneTreeTempName=geneTreeTempName);
  }else{
    distVecTable <- read.table(distVecName);
    distVec <- NULL
    for(i in 1:ncol(distVecTable)){
      distVec <- unlist(c(distVec, distVecTable[,i])); 
    }
  }
  if(saveXLVec){
    write(distVec, distVecName);
  }
  if(saveHistogram){
    hist(distVec,  xlab = 'expected extra lineages', ylab = "Frequency", main=NULL);
    dev.copy(pdf, paste("Histogram-",Exp.outputName, ".pdf", sep=""), height=height, width=width); #height=5.83, width=8.27 for A5 size
    dev.off();
  }
  if(only.XLvec == F){
    genetreeFiles <- list.files(path = genetree.Path, pattern = genetreePattern); # this will get the gene tree names
    genetreeFiles <- genetreeFiles[genetreeFiles != sptreeName];
    cat("A total of ", length(genetreeFiles), "gene trees were found:\n", genetreeFiles,"\n");
    myDist <- get.myGenetreeXL(genetreeFiles=genetreeFiles, genetree.Path=genetree.Path, sptree.dir=getwd(), sp.treeName=sptreeName, species=species, Imap=Imap.txt, 
                               msTaxaNames=F, taxaVec=NULL, PhyloImapName=PhyloImapName, PhylonetDir=PhylonetDir);
    likelihood <- getLikelihood(genetreeFiles = myDist, 
                                distVec=distVec, 
                                sp.treeName = sptreeName,
                                PhylonetDir = PhylonetDir, wd = wd, genetree.Path=genetree.Path,
                                Imap = Imap.txt, species = species, msTaxaNames=F, taxaVec = NULL, PhyloImapName=PhyloImapName);
    result <- data.frame(distVecName, gfDesc, likelihood);
    colnames(result) <- c("distVec.Source", "M", "likelihood");
    return(result);
  }
}

run.xl.pipeline.parallel <- function(wd=getwd(), model.table, cores=1, sptreeName="sp-tree.tre", genetreePattern="tre$",genetree.Path=getwd(), ImapName="Imap.txt",
                                     ms.exec="./ms", brlength.correction=NULL, nloci=10000, PhylonetDir=NULL,
                                     saveHistogram=T, height=8.27, width=8.27){
  library(foreach);
  library(doMC);
  library(ape);
  library(gmp);
  registerDoMC(cores);
  models <- unique(as.integer(model.table$model.number));
  foreach(i=1:length(models), .combine='rbind') %dopar%{
    subtable <- model.table[model.table$model.number == models[i],];
    time <- subtable$time;
    species1 <- as.character(subtable$species1);
    species2 <- as.character(subtable$species2);
    M <- subtable$M;
    PhyloImapName <- paste("PhyloImap", i, ".txt", sep="");
    geneTreeTempName <- paste("genetree", i, ".temp", sep="");
    run.xl.pipeline(wd=wd, sptreeName=sptreeName, genetreePattern=genetreePattern, genetree.Path=genetree.Path, ImapName=ImapName,
                    ms.exec=ms.exec, brlength.correction=brlength.correction, nloci=nloci, PhylonetDir=PhylonetDir,
                    PhyloImapName=PhyloImapName, geneTreeTempName=geneTreeTempName, saveXLVec=T, loadXLvec=F,
                    only.XLvec=T,
                    time=time, species1=species1, species2=species2, M=M,
                    saveHistogram=saveHistogram, height=height, width=width);
    
  }
  final.table <- data.frame();
  for(i in 1:length(models)){
    subtable <- model.table[model.table$model.number == models[i],];
    time <- subtable$time;
    species1 <- as.character(subtable$species1);
    species2 <- as.character(subtable$species2);
    M <- subtable$M;
    PhyloImapName <- paste("PhyloImap", i, ".txt", sep="");
    geneTreeTempName <- paste("genetree", i, ".temp", sep="");    
    table <- run.xl.pipeline(wd=wd, sptreeName=sptreeName, genetreePattern=genetreePattern, genetree.Path=genetree.Path, ImapName=ImapName,
                             ms.exec=ms.exec, brlength.correction=brlength.correction, nloci=nloci, PhylonetDir=PhylonetDir,
                             PhyloImapName=PhyloImapName, geneTreeTempName=geneTreeTempName, saveXLVec=F, 
                             loadXLvec=T,
                             only.XLvec=F,
                             time=time, species1=species1, species2=species2, M=M,
                             saveHistogram=F, height=height, width=width);
    
    final.table <- rbind(final.table, table);
  }
  return(final.table);
}
