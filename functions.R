library(ape);

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
  for(i in species){ # to replace species names by numbers, then ms can take them
    gf.code <- gsub(pattern=paste("([[:space:]])",i, "([[:space:]])", sep=""), # the space is just to be sure it's matching taxa names, even when taxa names are numbers.
                    replacement=paste("\\1",count,"\\2", sep=""), 
                    x=gf.code);
    count <- count+1;
  }
  gf.code <- gsub(pattern="(\\w+)[[:space:]]$", replacement="\\1", x=gf.code); # to delete space at the end of code
  gf.code <- paste("-em", gf.code);
  return(gf.code);
}

# get individual number from a fasta matrix and codify ms code
  # this function can be useful when having missing data, but this is not implemented in the tutorial provided
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

# execute ms
ms <- function(totalTaxa, nloci, species, taxaVec, ms.topology, gf.code = NULL, popSize.code = NULL, outputName){
  if(length(gf.code) >1){
    gf.code <- paste(gf.code, collapse=" ");
  }
  command <- paste("./ms", 
                   totalTaxa, 
                   nloci, 
                   "-T -I", 
                   length(species), 
                   taxaVec, 
                   topology,
                   gf.code,
                   popSize.code,
                   "| tail +4 | grep -v // >", 
                   outputName);
  system(command);
}

### write PhyloImap
writePhyloImap <- function(Imap, species, msTaxaNames = F, taxaVec = NULL){
  PhyloImap <- character();
    # writting Imap required by Phylonet
  if(msTaxaNames){
    if(length(taxaVec) == 0){
      stop("To write PhyloImap using ms taxa names, you need to provide a taxaVec\n");
    }
    taxaVec <- as.integer(unlist(strsplit(taxaVec, split=" ")));
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
  write(PhyloImap[1:length(species)], "PhyloImap.txt");
}

#### sp tree vs genetrees
sptree.vs.genetrees <- function(wd= getwd(), PhylonetDir = getwd(), 
                               tree.newick, genetreeNames, species, Imap, msTaxaNames = F, taxaVec = NULL, sp.treeName, sptree.dir = getwd()){
  writePhyloImap(Imap, species, msTaxaNames = msTaxaNames, taxaVec = taxaVec);
  distVec <- integer();
  cat("Reading expected gene trees from file:",genetreeNames,"\n");
  genetrees <- read.tree(genetreeNames);
  if(class(genetrees) == "multiPhylo"){
    count <- 1;
    for(l in 1:length(genetrees)){
      write.tree(genetrees[[l]], "genetree.temp");
      command <- paste("java -jar ", PhylonetDir, " deep_coal_count ", sptree.dir, "/", sp.treeName," ", wd, "/genetree.temp -a ", wd, "/PhyloImap.txt", sep="");
      deepcoal <- system(command, intern = T);
      deepcoal <- as.integer(gsub(pattern=".*lineages: (\\d+)", replacement="\\1", deepcoal[2]));
      distVec <- c(distVec, deepcoal);
      cat("Estimating expected distribution of extra lineages, progress: ",l , " / ", length(genetrees), "\n");
      count <- count +1;
    }
  }else{
    write.tree(genetrees, "genetree.temp");
    command <- paste("java -jar ", PhylonetDir, " deep_coal_count ", wd, "/", sp.treeName," ", wd, "/genetree.temp -a ", wd, "/PhyloImap.txt", sep="");
    deepcoal <- system(command, intern = T);
    deepcoal <- as.integer(gsub(pattern=".*lineages: (\\d+)", replacement="\\1", deepcoal[2]));
    distVec <- c(distVec, deepcoal);
  }
  return(distVec)
}
 
    
# comparing real gene trees with model
getLikelihood <- function(sptree, genetreeFile, distVec, sp.treeName = NULL,
                          dist.topo = T, 
                          Phylonet = F, PhylonetDir = getwd(), wd = getwd(),
                          Imap, species, taxaVec, sptree.dir=getwd()){
  library(gmp);
  pVec <- numeric();
  nVec <- integer();
  myDistVec <- integer();
  H <- length(distVec);
  N <- length(genetreeFiles);
  for(i in 1:N){
    mygenetreeFile <- file.path("genetrees",genetreeFiles[i]);
    writePhyloImap(Imap, species, msTaxaNames = F, taxaVec = NULL);
      # comparing with real data
    command <- paste("java -jar ", PhylonetDir, " deep_coal_count ", sptree.dir, "/", sp.treeName," ", wd, "/", mygenetreeFile, " -a PhyloImap.txt", sep="");
    mydeepcoal <- system(command, intern = T);
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
  my.freq <- unique(myDistVec);
  count <- 1;
  for(j in my.freq){
    nVec[count] <- sum(myDistVec == j);
    count <- count+1;
  }
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

### get.pvalue  - this function will return only the probability vector (could be needed when having missing data)
get.p <- function(sptree, genetreeFile, distVec, sp.treeName = NULL,
                          dist.topo = T, 
                          Phylonet = F, PhylonetDir = getwd(), dir = getwd(),
                          Imap, species, taxaVec){
  library(gmp);
  pVec <- numeric();
  nVec <- integer();
  myDistVec <- integer();
  H <- length(distVec);
  N <- length(genetreeFiles);
  for(i in 1:N){
    mygenetreeFile <- genetreeFiles[i];
    writePhyloImap(Imap, species, msTaxaNames = F, taxaVec = NULL);
    # comparing with real data
    command <- paste("java -jar ", PhylonetDir, " deep_coal_count ", dir, "/", sp.treeName," ", dir, "/", mygenetreeFile, " -a PhyloImap.txt", sep="");
    mydeepcoal <- system(command, intern = T);
    myDist <- as.integer(gsub(pattern=".*lineages: (\\d+)", replacement="\\1", mydeepcoal[2]));
    cat("Extra lineages observed:", myDist, "\n");
    if(all(myDist != myDistVec)){
      h <- sum(distVec == myDist);
      if(h == 0){ # this is to prevent p = 0, then we choose a very small number
        pVec <- c(pVec, 0.0000000001);
      }else{
        p <-  h/H;
        hist(distVec);
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
  # N!/ni! * pi^ni
  # ln N! - ln ni! + ni * ln pi
  calculation <-  sum(nVec*log(pVec)) - sum(log(factorialZ(nVec))); 
  return(calculation);
}

# this function gets the N vector
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

# get likelihood when providing a pVec, nVec and N
get.L <- function(nVec, pVec, N){
  # N!/ni! * pi^ni
  # ln N! - ln ni! + ni * ln pi
  likelihood <- log(factorialZ(N)) - sum(log(factorialZ(nVec))) + sum(nVec*log(pVec));
  return(likelihood);
}
