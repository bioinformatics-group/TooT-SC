firstrun=F
blastpSeq<- function(seq, start.pos = 1L, end.pos = nchar(seq), 
                     blastp.path = NULL, makeblastdb.path = NULL, 
                     database.path = NULL, silent = TRUE, 
                     evalue = 10L, output.path=normalizePath(resultspath)){
  
  if (Sys.which('makeblastdb') == '' & is.null(makeblastdb.path))
    stop('Please install makeblastdb (included in the NCBI BLAST+) or specify makeblastdb.path.')
  
  if (Sys.which('blastp') == '' & is.null(blastp.path))
    stop('Please install blastp (included in the NCBI BLAST+) or specify blastp.path.')
  
  makeblastdb.path = if (!is.null(makeblastdb.path)) makeblastdb.path else Sys.which('makeblastdb')
  blastp.path = if (!is.null(blastp.path)) blastp.path else Sys.which('blastp')
  
  if (is.null(database.path)) stop('Must specify the database (a FASTA file) path')
  if (is.null(output.path)) stop('Must specify the output path')
  
  N = end.pos - start.pos + 1L
  
  # Prepare data for Blastp
  cmddb = paste0(shQuote(makeblastdb.path), ' -dbtype prot -in ', 
                 shQuote(database.path),' -parse_seqids')        
  if(firstrun)
  {
    print("performing make blastdb")
    # print( cmddb)
    if (silent == TRUE) system(cmddb, ignore.stdout = TRUE) else system(cmddb)
  }
  
  # Basic parameters for Blastp
  tmp = tempfile('Blastp')
  queryFasta = paste0(tmp, '.fasta')
  tabularfile= paste0(tmp, '.txt')
  querySeq = Biostrings::AAStringSet(as.character(seq))
  Biostrings::writeXStringSet(querySeq, queryFasta)          
  # Additional parameters for Blastp
  if (!is.null(evalue)) {
    if (evalue <= 0L) {
      stop('evalue must be > 0')
    }
  }
  outFile <- file.path(output.path, "out.txt")
  # Run Blastp
  cmdblastp = paste(
    paste0(shQuote(blastp.path),
           ' -comp_based_stats 1 -db ', shQuote(database.path),
           ' -query ', shQuote(queryFasta),  ' -outfmt 6',' -out ', shQuote(outFile)))
  
  print("******************************")
  print(cmdblastp)
  if (silent == TRUE) system(cmdblastp, ignore.stdout = F) else system(cmdblastp)      
  #get the hit sequences Id
  if(file.info(outFile)$size != 0) # if there some hits are found
   {
    data = read.table(outFile)
    HomologousSeqIds= data$V2
    dindex <- which(duplicated(HomologousSeqIds))
    if(length(dindex) !=0 )
    {
      HomologousSeqIds=HomologousSeqIds[-dindex]# remove duplicates if any
    }
    # print(length(HomologousSeqIds))

    if(length(HomologousSeqIds)>=120)
      HomologousSeqIds<-HomologousSeqIds[1:120] 
    else
      HomologousSeqIds<- HomologousSeqIds

    fileName<-paste0(output.path,"H.txt")
    fileConn<-file(fileName)
    write(as.character(HomologousSeqIds), fileConn)
    close(fileConn)

    #get the cossponding Fasta file
    seqFile <- file.path(output.path, "seq.txt")
    getseqcmd= paste0(shQuote(Sys.which('blastdbcmd')),' -db ',shQuote(database.path), ' -entry_batch ', fileName, ' -out ', seqFile)
    # print(getseqcmd)
    if (silent == TRUE) system(getseqcmd, ignore.stdout = TRUE) else system(getseqcmd)
      
  }else{
    print (paste0( output.path,"---No hits found"))
  }
}

FilteredMSA= function(path)
{
  #print(path)
  setwd(path)
  file.remove("seq.tcs_column_filter4_fasta")
  
  tcsScorecmd<-paste0("t_coffee seq.txt -mode psicoffee -blast_server=LOCAL -protein_db ",dbpath2,"uniref50-tm.fasta  -output tcs_residue_filter3_fasta,clustalw_aln,tcs_column_filter4_fasta,score_html")
  system(tcsScorecmd)
  system('rm *.prf')
  #Removing columns with 80%gaps or more 
  removegapscmd= paste0("t_coffee -other_pg seq_reformat -in seq.tcs_column_filter4_fasta", " -output  fasta > filteredSeq.fasta")
  system(removegapscmd)
  
}
ClaculateCompositions_return<- function(subdirName, seqfile)
{ 
  
  print(seqfile)
  
  se<- read.fasta(seqfile,seqtype = "AA",as.string=T)
  
  seqlist<- unlist(lapply(se,as.character))
  seqlist<-sapply(seqlist,gsub,pattern="[^GPAVLIMCFYWHKRQNEDST]",replacement="") # getting rid of gaps? not the smartest thing to do, need more work
  names(seqlist)<- sub(".+\\:","",names(se))
  names(seqlist)<- sub("\\|.*","",sub(".+?\\|","", names(seqlist)))
  seqlist=seqlist[which(nchar(seqlist)>30)]

  MSADC<- lapply(seqlist, extractDC)
  outputDC <- matrix(unlist(MSADC), ncol = 400, byrow = TRUE)
  rownames(outputDC)<- names(seqlist)
  write.csv( outputDC,paste0(subdirName,"PAAC",".csv"))
  return(list(DC=outputDC))
  
}

PreparedataforMSAAAC = function(seq, start.pos = 1L, end.pos = nchar(seq), 
                                blastp.path = NULL, makeblastdb.path = NULL, 
                                database.path = NULL, iter = 5, silent = TRUE, 
                                evalue = 10L, output.path=normalizePath(resultspath)) {
  #1- run blastp on datafiles
  seqname=names(seq)
  SeqDirectory=file.path(output.path,names(seq))
  
  dir.create(SeqDirectory, showWarnings = FALSE, recursive = FALSE, mode = "0777")   
  blastpSeq(seq, start.pos , end.pos ,  blastp.path, makeblastdb.path , database.path , silent ,evalue, SeqDirectory)
  
  
  #2- Do Filtered MSA 
  
  #print("creating filtered Seq")
  
  #FilteredMSA(SeqDirectory)  
  
}
MSA_PAAC<- function(fastafile)
{
  seqs<- readFASTA(fastafile)
  names(seqs)<- sub("\\|.*","",sub(".+?\\|","", names(seqs)))
  for(j in c(1:length(seqs)))
  {
    x<- seqs[j]
    PreparedataforMSAAAC(seq= x,database.path=swissprotdb,output.path=intermediateFiles)
  }
  

AAfiles<-  names(seqs)
MSAPAAC <- matrix(ncol=400+1,nrow=length(AAfiles))
dfMSAPAAC <- data.frame(matrix(ncol = 400+1, nrow =0)) 

for(i in c(1:length(AAfiles)))
{
  subdirName <- file.path(intermediateFiles, AAfiles[i])
    seqFile <- file.path(subdirName, "seq.txt")
    if(file.exists(seqFile) && (file.info(seqFile)$size > 0))
    {
        AllComp<- ClaculateCompositions_return(subdirName, seqFile)
        outputDC <- AllComp$DC
        DC<- apply(outputDC, 2, mean)
        MSAPAAC[i,]<-   c("NA",DC)
      
    }else{ # if the sequence does not have any hits, just compute the PAAC from the single query
        print("seq.txt does not exist")
        DC<- extractDC(as.character(seqs[i]))
        MSAPAAC[i,]<-   c("NA",DC)
    }

}
dfMSAPAAC<- rbind(dfMSAPAAC,MSAPAAC)


write.csv(cbind(UniprotID=AAfiles,dfMSAPAAC), file = file.path(compositions,"MSA_PAAC.csv"),row.names = F)
}
