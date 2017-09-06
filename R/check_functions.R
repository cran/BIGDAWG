#' Check Input Parameters
#'
#' Check input parameters for invalid entries.
#' @param HLA Logical indicating whether data is HLA class I/II genotyping data only.
#' @param All.Pairwise Logical indicating whether all pairwise loci should be analyzed in haplotype analysis.
#' @param Trim Logical indicating if HLA alleles should be trimmed to a set resolution.
#' @param Res Numeric setting what desired resolution to trim HLA alleles.
#' @param EVS.rm Logical indicating if expression variant suffixes should be removed.
#' @param Missing Numeric setting allowable missing data for running analysis (may use "ignore").
#' @param Cores.Lim Interger setting the number of cores accessible to BIGDAWG (Windows limit is 1 core).
#' @param Return Logical Should analysis results be returned as list.
#' @param Output Logical Should analysis results be written to output directory.
#' @param Merge.Output Logical Should analysis results be merged into a single file for easy access.
#' @param Verbose Logical Should a summary of each analysis be displayed in console.
#' @note This function is for internal use only.
Check.Params <- function (HLA,All.Pairwise,Trim,Res,EVS.rm,Missing,Cores.Lim,Return,Output,Merge.Output,Verbose) {
  
  # Logicals: HLA=TRUE, All.Pairwise=FALSE, EVS.rm=FALSE, Trim=FALSE, Return=FALSE, Merge.Output=FALSE, Verbose=TRUE, Output=TRUE,
  # Numerics: Res=2, Missing=2, Cores.Lim=1L
  # Untested: Data, Results.Dir, Run.Tests, Loci.Set
  
  if( !is.logical(HLA) ) { Err.Log("P.Error","HLA") ; stop("Conversion Stopped.",call.=FALSE) }
  if( !is.logical(All.Pairwise) ) { Err.Log("P.Error","All.Pairwise") ; stop("Conversion Stopped.",call.=FALSE) }
  if( !is.logical(EVS.rm) ) { Err.Log("P.Error","EVS.rm") ; stop("Conversion Stopped.",call.=FALSE) }
  if( !is.logical(Trim) ) { Err.Log("P.Error","Trim") ; stop("Conversion Stopped.",call.=FALSE) }
  if( !is.logical(Return) ) { Err.Log("P.Error","Return") ; stop("Conversion Stopped.",call.=FALSE) }
  if( !is.logical(Merge.Output) ) { Err.Log("P.Error","Merge.Output") ; stop("Conversion Stopped.",call.=FALSE) }
  if( !is.logical(Verbose) ) { Err.Log("P.Error","Verbose") ; stop("Conversion Stopped.",call.=FALSE) }
  if( !is.logical(Output) ) { Err.Log("P.Error","Output") ; stop("Conversion Stopped.",call.=FALSE) }
  if( !is.numeric(Res) ) { Err.Log("P.Error","Res") ; stop("Conversion Stopped.",call.=FALSE) }
  if( !is.numeric(Missing) ) { Err.Log("P.Error","Missing") ; stop("Conversion Stopped.",call.=FALSE) }
  if( !is.numeric(Cores.Lim) || !is.integer(Cores.Lim) ) { Err.Log("P.Error","Cores.Lim") ; stop("Conversion Stopped.",call.=FALSE) }
  
}

#' Check Cores Parameters
#'
#' Check cores limitation for OS compatibility
#' @param Cores.Lim Integer How many cores can be used.
Check.Cores <- function(Cores.Lim) {
  if ( Cores.Lim!=1L ) {
    Cores.Max <- as.integer( floor( parallel::detectCores() * 0.9) )
    if(Sys.info()['sysname']=="Windows" && as.numeric(Cores.Lim)>1) {
      Err.Log("Windows.Cores") ; stop("Conversion stopped.",call. = F)
    } else if( Cores.Lim > Cores.Max ) { Cores <- Cores.Max
    } else { Cores <- Cores.Lim }
  } else { Cores <- Cores.Lim }
  return(Cores)
}

#' HLA Formatting Check
#' 
#' Checks data to see if HLA data is properly formatted.
#' @param x All columns of HLA genotyping data.
#' @note This function is for internal BIGDAWG use only.
CheckHLA <- function(x) {
  #Return TRUE if properly formatted HLA
  
  x[is.na(x)] <- "00:00" # NA cells
  x[x=="^"] <- "00:00" # absent cells
  test <- apply(x,MARGIN=c(1,2),FUN=function(z) length(unlist(strsplit(as.character(z),split=":"))))
  test <- apply(test,MARGIN=2,FUN=min)
  Flag <- as.logical(min(test)==2)
  
  return(Flag)
}

#' Allele Name Format Fix
#' 
#' Separate locus and allele names if data is formatted as Loci*Allele
#' @param Output Logical indicating if Error logging should be written to a file.
#' @param Tab All columns of HLA genotyping data.
#' @note This function is for internal BIGDAWG use only. 
FixAlleleName <- function(Output,Tab) {
  fixCell <- apply(Tab[,3:ncol(Tab)],MARGIN=c(1,2),FUN=function(x) grepl("\\*",na.omit(x)))
  if( sum(fixCell)>0 ) {
    Tab[Tab=="^"] <- "*^"
    Format.Check <- length(apply(Tab[,3:ncol(Tab)],MARGIN=2,FUN=function(x) which(grepl("\\*",na.omit(x))==F)))
    if( Format.Check>0 ) {
      Err.Log(Output,"Uneven.Prefix")
      stop("Analysis Stopped.",call. = F)
    }
    Tab[,3:ncol(Tab)] <- apply(Tab[,3:ncol(Tab)],MARGIN=c(1,2),FUN=function(x) unlist(lapply(strsplit(x,split="\\*"),"[",2)))
  }
  return(Tab)
}

#' Loci Legitimacy Check
#' 
#' Checks available loci against data to ensure complete overlap.
#' @param x Loci available in exon protein list alignment object.
#' @param y Unique column names
#' @note This function is for internal BIGDAWG use only.
CheckLoci <- function(x,y) {
  #Returns TRUE if absent locus(loci) encountered
  #x=Loci available in ExonPtnList
  #y=Loci.Set from data
  
  Output <- list()
  
  y <- unique(unlist(y))
  
  Flag <- ( !sum(y %in% x) == length(y) )
  Output[['Flag']] <- Flag
  if(Flag) {
    Output[['Loci']] <- paste(y[!y %in% x],collapse=",")
  } else {
    Output[['Loci']] <- NA
  }
  
  return(Output)
}

#' Allele Legitimacy Check
#' 
#' Checks available alleles against data to ensure complete overlap.
#' @param x Exon protein list alignment object.
#' @param y Genotypes from data file
#' @param z1 loci in data file
#' @param z2 Genotype column names
#' @note This function is for internal BIGDAWG use only.
CheckAlleles <- function(x,y,z1,z2) {
  #Returns TRUE if unknown allele(s) encountered
  #Checks at 3 levels of resolution: Full, 2-Field, and 1-Field
  #x=ExonPtnList
  #y=genos
  #z1=loci
  #z2=loci.ColNames
  
  Output <- list()
  for(i in 1:length(z1)) {
    
    Locus.tmp <- z1[i]
    
    # Available Alleles
    x.locus <- x[[Locus.tmp]]
    
    # Data Alleles
    y.locus <- y[which(z2==Locus.tmp)]
    y.locus <- matrix(c(unlist(y.locus[,1]),y.locus[,2]),ncol=1)
    y.locus <- na.omit(y.locus)
    
    # Remove Absent Allele Calls
    Allele.rm <- which(y.locus=="^")
    if(length(Allele.rm)>0) { y.locus <- y.locus[-Allele.rm] }
    
    #identify minimum resolution
    tmp <- sapply(y.locus,FUN=strsplit,split=":")
    tmp <- lapply(tmp,unlist)
    Res <- min(unlist(lapply(tmp,length)))
    
    y.alleles <- sort(unique(unlist(y.locus)))
    x.alleles <- unique(x.locus[,'Allele'])
    
    #format according to minimum resolution
    if (Res>=2) {
      y.alleles <- unique(sapply(y.alleles,GetField,Res=2))
      x.alleles <- unique(x.locus[,'Trimmed'])
    } else if (Res==1) {
      y.alleles <- unique(sapply(y.alleles,GetField,Res=1))
      x.alleles <- unique(sapply(unique(x.locus[,'Trimmed']),GetField,Res=1))
    }
    
    Output[[Locus.tmp]] <- list(Flag=ifelse(!sum(y.alleles %in% x.alleles)==length(y.alleles),T,F),
                                Allele=paste(Locus.tmp,paste(y.alleles[!y.alleles %in% x.alleles],collapse=","),sep="*"))
    
    i = i + 1
    Output
    
  }
  return(Output)
}

#' Data Summary Function
#' 
#' Summary function for sample population within data file.
#' @param Tab Loci available in exon protein list alignment object.
#' @param All.ColNames Column names from genotype data.
#' @param rescall HLA resolution set for analysis.
#' @param HLA HLA bigdawg argument passed to function
#' @param Verbose Summary display carryover from BIGDAWG function.
#' @param Output Data output carryover form BIGDAWG function
#' @note This function is for internal BIGDAWG use only.
PreCheck <- function(Tab,All.ColNames,rescall,HLA,Verbose,Output) {
  
  Grp0 <- which(Tab[,2]==0)
  Grp1 <- which(Tab[,2]==1)
  nGrp0 <- length(Tab[Grp0,2])
  nGrp1 <- length(Tab[Grp1,2])
  
  if(min(nGrp0,nGrp1)==0) {
    Err.Log(Output,"Case.Con")
    stop("Analysis Stopped.",call. = F)
  }
  
  Loci <- as.list(unique(All.ColNames[3:length(All.ColNames)]))
  nLoci <- length(Loci)
  GTYPE <- Tab[,3:ncol(Tab)]
  colnames(GTYPE) <- All.ColNames[3:length(All.ColNames)]
  nGTYPE <- unlist(lapply(Loci,function(x) length(unique(unlist(GTYPE[,which(colnames(GTYPE)==x)])))))
  Grp0un <- unlist(lapply(Loci,function(x) length(unique(unlist(GTYPE[Grp0,which(colnames(GTYPE)==x)])))))
  Grp1un <- unlist(lapply(Loci,function(x) length(unique(unlist(GTYPE[Grp1,which(colnames(GTYPE)==x)])))))
  
  nMissing <- unlist(lapply(Loci,function(x) sum(is.na(GTYPE[,which(colnames(GTYPE)==x)]))))
  Grp0miss <- unlist(lapply(Loci,function(x) sum(is.na(GTYPE[Grp0,which(colnames(GTYPE)==x)]))))
  Grp1miss <- unlist(lapply(Loci,function(x) sum(is.na(GTYPE[Grp1,which(colnames(GTYPE)==x)]))))
  
  if(Verbose) {
    cat("  Sample Summary\n")
    cat("    Sample Size (n):",nrow(Tab),"\n")
    cat("    ...Number of Controls/Cases:",paste(paste(nGrp0,nGrp1,sep="/"),collapse=", "),"\n")
    cat("    Allele Count (2n):",nrow(Tab)*2,"\n")
    cat("    Total loci in file:",nLoci,"\n")
    cat("    Unique loci:",paste(Loci,collapse=", "),"\n") 
    cat("    Unique alleles per locus:",paste(nGTYPE,collapse=", "),"\n")
    cat("    ...Unique in Controls/Cases:",paste(paste(Grp0un,Grp1un,sep="/"),collapse=", "),"\n")
    cat("    Missing alleles per locus:",paste(nMissing,collapse=", "),"\n")
    cat("    ...Missing in Controls/Cases:",paste(paste(Grp0miss,Grp1miss,sep="/"),collapse=", "),"\n")
    cat("\n")
  }
  
  if(HLA) {
    
    Grp0res <- max(unlist(lapply(Loci,function(x) max(unlist(lapply(strsplit(unlist(GTYPE[Grp0,which(colnames(GTYPE)==x)]),split=":"),length))))))
    Grp1res <- max(unlist(lapply(Loci,function(x) max(unlist(lapply(strsplit(unlist(GTYPE[Grp1,which(colnames(GTYPE)==x)]),split=":"),length))))))
    
    if(max(Grp0res,Grp1res)>4) {
      Err.Log(Output,"High.Res")
      stop("Analysis Stopped.",call. = F)
    }
    
    if(Verbose){
      cat("  Observed Allele Resolution\n")
      cat("    Max Resolution Controls:",paste(Grp0res,"-Field",sep=""),"\n")
      cat("    Max Resolution Cases:",paste(Grp1res,"-Field",sep=""),"\n")
      cat("    Defined Resolution:",rescall,"\n")
      if(Grp0res!=Grp1res){ cat("  ***** Warning. \n") }
      if(Grp0res!=Grp1res){ cat("  ***** There may exist a Case-Control field resolution imbalance.\n") }
      if(Grp0res!=Grp1res){ cat("  ***** Considering trimming to",paste(min(Grp0res,Grp1res),"-Field resolution.",sep=""),"\n") }
      cat("\n")
    }
  }
  
  if(HLA) {
    Out <- list(SampleSize=nrow(Tab),
                   No.Controls=nGrp0,
                   No.Cases=nGrp1,
                   AlleleCount=nrow(Tab)*2,
                   TotalLoci=nLoci,
                   Loci=paste(Loci,collapse=", "),
                   AllelePerLocus=paste(nGTYPE,collapse=", "),
                   MissingPerLocus=paste(nMissing,collapse=", "),
                   MaxResGrp0=paste(Grp0res,"-Field",sep=""),
                   MaxResGrp1=paste(Grp1res,"-Field",sep=""),
                   SuggestedRes=paste(min(Grp0res,Grp1res),"-Field",sep=""),
                   SetRes=rescall)
  } else {
    Out <- list(SampleSize=nrow(Tab),
                   No.Controls=nGrp0,
                   No.Cases=nGrp1,
                   AlleleCount=nrow(Tab)*2,
                   TotalLoci=nLoci,
                   Loci=paste(Loci,collapse=", "),
                   AllelePerLocus=paste(nGTYPE,collapse=", "),
                   MissingPerLocus=paste(nMissing,collapse=", "),
                   SetRes=rescall)
  }
  
  return(do.call(rbind,Out))
  
}
