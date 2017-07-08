#' Update function for protein aligment upon new IMGT HLA data release
#'
#' This updates the protein aligment used in checking HLA loci and alleles as well as in the amino acid analysis. Alignment must exist in database (ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/) or update will fail.
#' @param Add.Loci Character string or vector of loci that should be added to default loci (default = HLA-A,B,C,DRB1/3/4/5,DQA1,DQB1,DPA1,DPB1).
#' @param Restore Logical specifying if the original alignment file be restored.
#' @param Force Logical specifiying if update should be forced.
#' @param Output Logical indicating if error reporting should be written to file.
UpdateRelease <- function(Add.Loci=NULL,Force=F,Restore=F,Output=F) {
  
  if( !inherits(try(XML::readHTMLTable("http://cran.r-project.org/web/packages/BIGDAWG/index.html",header=F),silent=T),"try-error") ) {
  
    MainDir <- getwd()
    on.exit(setwd(MainDir), add = TRUE)
  
    getDir <- path.package('BIGDAWG')
    putDir <- paste(getDir,"/data",sep="")
    if(!dir.exists(putDir)) { dir.create(putDir) }
    
    if(!Restore) {
      
      #Check current version against BIGDAWG version
      if(!Force) {
        RV <- XML::readHTMLTable("http://www.ebi.ac.uk/ipd/imgt/hla/docs/release.html",header=T)
        RV.current <- as.character(lapply(RV,"[",1)[[1]][1,])
        RV.BIGDAWG <- unlist(strsplit(as.character(ExonPtnList$Release[[1]]),":"))[2]
        cat("Versions:\n","IMGT/HLA current: ",RV.current,"\n BIGDAWG version: ",RV.BIGDAWG,"\n")
        if(grepl(RV.current,RV.BIGDAWG)) { Flag <- T } else { Flag <- F }
      } else {
        Flag <- F
      }# End if() for setting Flag
      
      #Run Update if Flag = T
      if(Flag) {
  
        cat("\nYour database seems up to date. Use Force = T to force the update.")
        
      } else {
        
        # For creating UpdatePtnAlign.RData object 
        # Define download directory
          setwd(putDir)
          Safe <- dir()
          Safe <- c(Safe[!grepl(".txt",Safe)],"UpdatePtnAlign.RData")
        
        #STEP 1: Define Loci and Read in Reference Exon Map Files
          # Loci
          if(is.null(Add.Loci)) {
            Loci <- c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1","DRB3","DRB4","DRB5")
          } else {
            Loci <- unique(c(Add.Loci,"A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1","DRB3","DRB4","DRB5"))
          }
        
          # Map
          RefTab <- BIGDAWG::ExonPtnList$RefExons
          
          #########################################################################
          # Need to remove following lines if DRB split into single locus files
          # Currently DRB1/3/4/5 are contained in single DRB.prot file.
          Loci.rep <- Loci[-grep("DRB",Loci)]
          Loci.rep <- sort(c(Loci.rep,"DRB"))
          #########################################################################
          
        #STEP 2: Download protein alignments and other ancillary files
          cat("Updating reference object for the amino acid analysis.\n")
          cat("Downloading alignment files from the IMGT/HLA.\n")
          GetFiles(Loci.rep)
          Release <- read.table('Release.txt',sep="\t") # created during GetFiles download
        
        #STEP 3: Format alignments for exons of interest
          cat("Formatting alignment files.\n")
          for(i in 1:length(Loci)) { Locus <- Loci[i] ; ExonPtnAlign.Create(Locus,RefTab) }
        
        #STEP 4: Create ExonPtnAlign list object for BIGDAWG package
          AlignObj.Update(Loci,Release,RefTab)
        
        #STEP 5: Clean up
          cat("Cleaning up.\n")
          invisible(file.remove(dir()[which(dir() %in% Safe!=T)]))
        
          cat("Updated.\n")
      }
      
    } else if (Restore) {
      
      setwd(putDir)
      if(!file.exists('UpdatePtnAlign.RData')) { stop("No prior update to restore.", call.= F) }
      cat("Restoring original alignment reference object for amino acid analysis.\n")
      invisible(file.remove('UpdatePtnAlign.RData'))
      cat("Restored.\n")
      
    }
    
  } else {
    
    Err.Log(Output,"No.Internet")
    stop("Analysis stopped.",call.=F)
    
  }

}
