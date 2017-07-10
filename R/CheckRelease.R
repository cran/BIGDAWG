#' Function to Check Release Versions
#'
#' This updates the protein aligment used in checking HLA loci and alleles as well as in the amino acid analysis.
#' @param Package Logical to check for BIGDAWG package versions
#' @param Alignment Logical to check the IMGT/HLA database version for the alignment bundled with BIGDAWG.
#' @param Output Should any error be written to a file
#' @note Requires active internet connection.
CheckRelease <- function(Package=T,Alignment=T,Output=F) {
  
  if( !inherits(try(XML::readHTMLTable("http://cran.r-project.org/web/packages/BIGDAWG/index.html",header=F),silent=T),"try-error") ) {
  
    if(Package) {
      
      CranR <- as.character(XML::readHTMLTable("http://cran.r-project.org/web/packages/BIGDAWG/index.html",header=F)[[1]][1,2])
      GitHubR <- read.table("https://raw.githubusercontent.com/pappasd/BIGDAWG/master/NEWS",sep="\t",stringsAsFactors=F,nrows=1)
      GitHubR <- gsub("v","",unlist(strsplit(as.character(GitHubR),split=" "))[2])
      CurrR <- as.character(packageVersion('BIGDAWG') )
      
    }
    
    if(Alignment) {
      
      RV <- XML::readHTMLTable("http://www.ebi.ac.uk/ipd/imgt/hla/docs/release.html",header=T)
      RV.current <- as.character(lapply(RV,"[",1)[[1]][1,])
      
      UPL <- paste(path.package('BIGDAWG'),"/data/UpdatePtnAlign.RData",sep="")
      UpdatePtnList <- NULL ; rm(UpdatePtnList)
      if( file.exists(UPL) ) { 
        load(UPL)
        EPL <- UpdatePtnList
        rm(UpdatePtnList,UPL)
        UPL.flag=T
      } else { 
        EPL <- ExonPtnList
        UPL.flag=F }
      
      RV.BIGDAWG <- gsub(" ","",unlist(strsplit(as.character(EPL$Release[[1]]),":"))[2])
      
    }
   
    cat("\n")
    if(Package) { cat("BIGDAWG Versions:\n","Installed Version: ",CurrR,"\n CRAN Release Version: ",CranR,"\n Developmental version: ",GitHubR,"\n") }
    if(Package & Alignment) { cat("\n") }
    if(Alignment) { 
      if(UPL.flag) {
          cat("IMGT/HLA Versions:\n","IMGT/HLA Version: ",RV.current,"\n Installed version (from update): ",RV.BIGDAWG,"\n") 
        } else {
          cat("IMGT/HLA Versions:\n","IMGT/HLA Version: ",RV.current,"\n Installed version: ",RV.BIGDAWG,"\n") 
        }
      }
    cat("\n")
  
  } else {
    
    Err.Log(Output,"No.Internet")
    stop("Analysis stopped.",call.=F)
    
  }
    
}