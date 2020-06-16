#! /usr/bin/Rscript

###################################################
## Name: TooT_SC_V1_11.R
## input: fasta file containing the unknown/testing protein sequences
## output: The predicted class of each protein sequence, and the classes probabilities in csv format
## Author: Munira Alballa
##################################################

standardized<-function(x,rmean,rsd){((x-rmean)/rsd)}
pop.sd<-function(x){sqrt(sum((x-mean(x))^2)/length(x))}
normalize<- function(matrix){
  data=matrix
  standardizedData<- matrix
  for( i in 1:length(data[,1]) )#until L
  {
    #compute mean across the 20 aa
    rmean= mean(data[i,])
    rsd=pop.sd(data[i,])
    standardizedData[i,]<- standardized(data[i,],rmean,rsd)
  }
  return(standardizedData)
  
  
}
oneVsRestSVM<- function(testdata){
  probabilities<- matrix(data=0, nrow=length(testdata[,1]), ncol=length(substates))
  colnames(probabilities)<-as.numeric(factor(substates, levels = substates))
  for( z in 1:length(substates))
  {
    svm.fit<- readRDS(paste0(TooTSCdir,"/models/MSAPAAC_class",z,"_1vsall.rds"))
    
    svm.predtest<-predict(svm.fit,standardizedData, probability=T)
    probabilities[,z]<- attr(svm.predtest,"probabilities")[,"1"]
    
  }
  
  
  pred<-as.vector(apply(probabilities, 1, function(x) names(which(x == max(x)))[1]))
  return(cbind.data.frame(pred=pred,probabilities ,stringsAsFactors=F ))
  
}




args <- commandArgs(trailingOnly=TRUE)

terminate <- FALSE

out <- "."
TooTSCdir <- "."
db<-"./db/"
for(i in args){
  arg = strsplit(i, "=")[[1]];
  
  switch(arg[1],
         "-query"={
           query <- arg[2]
         },
         "-out"={
           out <- arg[2]
         },
         "-TooTSC"={
           TooTSCdir <- normalizePath(arg[2])
         },
         "-db"={
           db <- normalizePath(arg[2])
         },
         "-help"={
           cat("TooTSC v1.0 (Oct. 2019)\n")
           cat("\n")
           cat("Usage: TooTSC -query=<input> [-TooTSC=<TooTSCdir>] [-out=<outdir>] [-db=<database path>]\n")
           cat("\n")
           cat("\t<input> is your sequence input file in fasta format\n")
           cat("\t<out> is the output directory where you want the predicted results, formatted as csv\n")
           cat("\t\t<out> defaults to '",out,"'\n")
           cat("\t<TooTSCdir> is the directory where the base TooT-SC files are located")
           cat("\t\t<TooTSCdir> defaults to '",TooTSCdir,"'\n")
           cat("\t <database path> is the relative path to the database\n")
           cat("\t\t<database path> defaults to",TooTSCdir,"/db/\n")
           cat("\n")
           terminate <- TRUE
           break
         }
  )
}

if(!terminate) {
  
  if(!exists("query")) {
    stop("-query has not been passed")
  }



  test_fasta <- normalizePath(path.expand(query))
  resultspath <- paste0(normalizePath(path.expand(out)),"/")
  

  require(seqinr)
  library("Biostrings")
  library("stringr")
  require(protr)
  library(ISLR)
  library(e1071)
  library(caret)
  library(R.utils)
  wd=normalizePath(path.expand(".")) # change the the tool directory
  
  if (isAbsolutePath(db)){
    dbpath <- db
  }else{
    dbpath <- paste0(TooTSCdir, db)
  }
  
  #dbpath=paste0(TooTSCdir, "/db/")
  compostions=paste0(TooTSCdir,"/intermediate_files/Compositions/")
  intermediateFiles=paste0(TooTSCdir,"/intermediate_files/")
  substates<- c("Nonselective",
                "water",
                "inorganic cation",
                "inorganic anion",
                "organic anion",
                "organooxogyn",
                "amino acid and derivatives",
                "other Organonitrogen compound",
                "nucleotide",
                "Organic heterocyclic",
                "Miscellaneous" )
  
  
  #testing data with unknown substrates
  source(paste0(TooTSCdir,"/src/MSA_PAAC_git.R"))
  
  MSA_PAAC(test_fasta)
  testfeatuers = read.csv(paste0(compostions,"MSA_PAAC.csv"),sep=",")
  
  
  #normalize
  standardizedData <- normalize(as.matrix(testfeatuers[,c(-1,-2)]))
  #predict
  svmpred<- oneVsRestSVM(standardizedData)
  
  # write results
  seqs<- readFASTA(test_fasta)
  names(seqs)<- sub("\\|.*","",sub(".+?\\|","", names(seqs)))
  print(paste0( "Toot-SC output is found at: ", resultspath, "TooTSCout.csv"))
  write.csv(cbind(UniProtID=names(seqs),svmpred ),paste0(resultspath,"TooTSCout.csv"))
  
  
}
