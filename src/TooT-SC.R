#! /usr/bin/Rscript

suppressMessages(suppressWarnings(library(seqinr)))
suppressMessages(suppressWarnings(library("Biostrings")))
suppressMessages(suppressWarnings(library("stringr")))
suppressMessages(suppressWarnings(library(protr)))
suppressMessages(suppressWarnings(library(ISLR)))
suppressMessages(suppressWarnings(library(e1071)))
suppressMessages(suppressWarnings(library(caret)))
suppressMessages(suppressWarnings(library(R.utils)))
suppressMessages(suppressWarnings(library("funr")))

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
  probabilities<- matrix(data=0, nrow=length(testdata[,1]), ncol=length(substrates))
  colnames(probabilities)<-as.numeric(factor(substrates, levels = substrates))
  for( z in 1:length(substrates))
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

out <- normalizePath(".")
TooTSCdir <- normalizePath(file.path(funr::get_script_path(), ".."))
db <- normalizePath(file.path(TooTSCdir, "db"), mustWork = FALSE)
work <- normalizePath(".")
for(i in args){
  arg = strsplit(i, "=")[[1]];
  
  switch(arg[1],
         "-query"={
           tootscquery <- normalizePath(arg[2])
         },
         "-out"={
           out <- normalizePath(arg[2])
         },
         "-TooTSC"={
           TooTSCdir <- normalizePath(arg[2])
         },
         "-db"={
           db <- normalizePath(arg[2])
         },
         "-work"={
           work <- normalizePath(arg[2])
         },
         "-help"={
           cat("TooTSC v1.1 (Apr. 2022)\n")
           cat("\n")
           cat("Usage: TooTSC -query=<input> [-out=<outdir>] [-db=<database path>] [-work=<work path>] [-TooTSC=<TooTSCdir>]\n")
           cat("\n")
           cat("\t<input> is your sequence input file in fasta format\n")
           cat("\t<out> is the output directory where you want the predicted results, formatted as csv\n")
           cat("\t\t<out> defaults to '.' ('",out,"')\n", sep="")
           cat("\t<database path> is the path to the database\n")
           cat("\t\t<database path> defaults to '",db,"'\n", sep="")
           cat("\t<work path> is the path to the working directory for intermediate files. It will be created as needed.\n")
           cat("\t\t<database path> defaults to '.' ('",work,"')\n", sep="")
           cat("\t<TooTSCdir> is the directory where the base TooT-SC files are located\n")
           cat("\t\t<TooTSCdir> defaults to '",TooTSCdir,"'\n", sep="")
           cat("\n")
           terminate <- TRUE
           break
         }
  )
}

if(!terminate) {
  

#
# Validate that the query exists
#
  if(!exists("tootscquery")) {
    stop("-query has not been passed")
  }

if(!file.exists(tootscquery)) {
   stop("The specified query file does not exist: '", tootscquery,"'", sep="")
}

#
# Validate that the db directory and required db files exist
#
if(!file.exists(db)) {
   stop("The specified database directory does not exist: '", db,"'", sep="")
}


dbFiles <- c("SwissOct18.fasta.psi", "SwissOct18.fasta.psd", "SwissOct18.fasta.pog", "SwissOct18.fasta.psq", "SwissOct18.fasta.pin", "SwissOct18.fasta", "SwissOct18.fasta.phr")
missingFiles <- list()
for(file in dbFiles) {
   if(!file.exists(file.path(db, file))) {
      missingFiles <- append(missingFiles, file)
   }
}

if(length(missingFiles) > 0) {
   stop("Unable to find some files in your db directory ('", db,"')\n", paste(missingFiles, collapse=", "))
}

swissprotdb <- file.path(db, "SwissOct18.fasta");


#
# Validate the outpit dir
#
if(!file.exists(out)) {
   stop("The specified output directory does not exist: '", out,"'", sep="")
}


#
# Validate and set up the working directory
#
if(!file.exists(work)) {
   stop("The specified base for your working directory does not exist: '", work,"'", sep="")
}

if(!file.exists(file.path(work, "work"))) {
   dir.create(file.path(work, "work"))
}

intermediateFiles = file.path(work, "work", "TooT-SC")
if(!file.exists(intermediateFiles)) {
   dir.create(intermediateFiles)
}

compositions = file.path(intermediateFiles, "Compositions")
if(!file.exists(compositions)) {
   dir.create(compositions)
}

#
# Lastly, validate the TooTSC Dir that they might have passed... I hate this param. We should also validate that all other required pieces are there
# which should catch broken installs
#
if(!file.exists(TooTSCdir)) {
   stop("The specified base for your application does not exist does not exist (you may want to consider using the dynamically generated one): '", TooTSCdir,"'", sep="")
}

MSAPAACSource <- file.path(TooTSCdir, "src", "MSA_PAAC_git.R")
if(!file.exists(MSAPAACSource)) {
   stop("A required source file to run TooTSC does not exist, though it should be located next to this script: '", MSAPAACSource,"'", sep="")
}



  substrates<- c("Nonselective",
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

missingModels <- list()
  for( z in 1:length(substrates)) {
    modelFile <- paste0(TooTSCdir,"/models/MSAPAAC_class",z,"_1vsall.rds")
    if(!file.exists(modelFile)) {
      missingModels <- append(missingModels, modelFile)
    }
  }
if(length(missingModels) > 0) {
   stop("Unable to find some models in your models directory ('", file.path(TooTSCdir, "models"),"')\n", paste(missingModels, collapse=", "))
}

  #testing data with unknown substrates
  source(MSAPAACSource)
  
  MSA_PAAC(tootscquery)
  testfeatuers = read.csv(file.path(compositions,"MSA_PAAC.csv"),sep=",")
  
  
  #normalize
  standardizedData <- normalize(as.matrix(testfeatuers[,c(-1,-2)]))
  #predict
  svmpred<- oneVsRestSVM(standardizedData)
  svmpred$pred<- substrates[as.numeric(svmpred$pred)]
  names(svmpred)<-c("pred",paste(substrates,"probability") )
  
  # write results
  seqs<- readFASTA(tootscquery)
  names(seqs)<- sub("\\|.*","",sub(".+?\\|","", names(seqs)))
  print(paste0( "Toot-SC output is found at: ", file.path(out,"TooTSCout.csv")))
  write.csv(cbind(UniProtID=names(seqs),svmpred ), file.path(out,"TooTSCout.csv"))
  
}
