#' BIGDAWG wrapper function
#'
#' This is the main wrapper function for each analysis.
#' @param Data Name of the genotype data file.
#' @param HLA Logical indicating whether data is HLA class I/II genotyping data only.
#' @param Run.Tests Specifics which tests to run.
#' @param Loci.Set Input list defining which loci to use for analyses (combinations permitted). 
#' @param All.Pairwise Logical indicating whether all pairwise loci should be analyzed in haplotype analysis.
#' @param Trim Logical indicating if HLA alleles should be trimmed to a set resolution.
#' @param Res Numeric setting what desired resolution to trim HLA alleles.
#' @param EVS.rm Logical indicating if expression variant suffixes should be removed.
#' @param Missing Numeric setting allowable missing data for running analysis (may use "ignore").
#' @param Cores.Lim Interger setting the number of cores accessible to BIGDAWG (Windows limit is 1 core).
#' @param Results.Dir Optional, string of full path directory name for BIGDAWG output.
#' @param Return Logical Should analysis results be returned as list.
#' @param Output Logical Should analysis results be written to output directory.
#' @param Merge.Output Logical Should analysis results be merged into a single file for easy access.
#' @param Verbose Logical Should a summary of each analysis be displayed in console.
#' @examples
#' ### The following examples use the synthetic data set bundled with BIGDAWG
#' 
#' # Haplotype analysis with no missing genotypes for two loci sets
#' # Significant haplotypic association with phenotype
#' # BIGDAWG(Data="HLA_data", Run.Tests="H", Missing=0, Loci.Set=list(c("DRB1","DQB1")))
#' 
#' # Hardy-Weinberg and Locus analysis ignoring missing data
#' # Significant locus associations with phenotype at all but DQB1
#' # BIGDAWG(Data="HLA_data", Run.Tests="L", Missing="ignore")
#' 
#' # Hardy-Weinberg analysis trimming data to 2-Field resolution
#' # Significant locus deviation at DQB1
#' BIGDAWG(Data="HLA_data", Run.Tests="HWE", Trim=TRUE, Res=2)
BIGDAWG <- function(Data, HLA=TRUE, Run.Tests, Loci.Set, All.Pairwise=FALSE, Trim=FALSE, Res=2, EVS.rm=FALSE, Missing=2, Cores.Lim=1L, Results.Dir, Return=FALSE, Output=TRUE, Merge.Output=FALSE, Verbose=TRUE) {

  options(warn=-1)
  
  MainDir <- getwd()
  on.exit(setwd(MainDir), add = TRUE)
    
  cat(rep("=",40))
  cat("\n          BIGDAWG: Bridging ImmunoGenomic Data Analysis Workflow Gaps\n")
  cat(rep("=",40),"\n")
  
  # Version Checking
  # BIGDAWG::CheckRelease(Alignment=F,Output=F)
  # if(HLA==T) { BIGDAWG::CheckRelease(Package=F,Output=F) }
  
  #Define Output Directory
  if (Output) {
    if(missing(Results.Dir)) {
      OutDir <- paste(MainDir,"/output ",format(Sys.time(), "%d%m%y %H%M%S"),sep="")
      dir.create(OutDir)
    } else {
      OutDir <- Results.Dir
    }
  }
  
  #Check for the updated ExonPtnList 'UpdatePtnList' and use if found.
  UpdatePtnList <- NULL ; rm(UpdatePtnList)
  UPL <- paste(path.package('BIGDAWG'),"/data/UpdatePtnAlign.RData",sep="")
  if( file.exists(UPL) ) { 
    load(UPL)
    EPL <- UpdatePtnList
    rm(UpdatePtnList)
    UPL.flag=T
  } else { 
    EPL <- BIGDAWG::ExonPtnList
    UPL.flag=F }
  
  cat("\n>>>>>>>>>>>>>>>>>>>>>>>>> BEGIN Analysis <<<<<<<<<<<<<<<<<<<<<<<<<\n\n")
  
  # Output object
  BD.out <- list()
  
# ===================================================================================================================================== ####
# Read in Data ________________________________________________________________________________________________________________________ ####
  
  NAstrings=c("NA","","****","-","na","Na")
  
  if (Data=="HLA_data") {
    Tab <- BIGDAWG::HLA_data
    colnames(Tab) <- toupper(colnames(Tab))
    All.ColNames <- gsub(".1","",colnames(Tab),fixed=T)
    rownames(Tab) <- NULL
    DRBFLAG <- NULL
    if(Output) { setwd(OutDir) }
    
  } else {
    
    # Read in data and Pre-process
    if(!file.exists(Data)) {
      Err.Log(Output,"Bad.Filename", Data)
      stop("Analysis stopped.",call.=F) }
    Tab <- read.table(Data, header = T, sep="\t", stringsAsFactors = F, na.strings=NAstrings, fill=T, comment.char = "#", strip.white=T, blank.lines.skip=T, colClasses="character")
    if(HLA==T) { if(sum(grepl("DRB3.4.5",colnames(Tab)))>0) { colnames(Tab) <- gsub("DRB3.4.5","DRB345",colnames(Tab))  } }
    colnames(Tab) <- toupper(colnames(Tab))
    All.ColNames <- gsub(".1","",colnames(Tab),fixed=T)
    rownames(Tab) <- NULL
    Tab <- rmABstrings(Tab)
    
    if(Output) { setwd(OutDir) }
    
    # Separate DRB345 if exists as single column pair and check zygosity
    if(HLA==T) {
      if(sum(grepl("DRB345",colnames(Tab)))>0) {
        
        cat("Processing DRB3/4/5 column data.\n")
        
        DRBFLAG <- T
        getCol <- grep("DRB345",colnames(Tab))
        
        # Stop if not Locus*Allele formatting
        if( sum(grepl("\\*",Tab[,getCol]))==0 ) {
          Err.Log(Output,"Bad.DRB345.format")
          stop("Analysis Stopped.",call. = F)
        }
        
        # Expand DRB3/4/5 to separate column pairs
        Tab <- DRB345.parser(Tab)
        
        # Check for locus hemizygosity by DR haplotype
        Tab.list <- lapply(seq_len(nrow(Tab)),FUN=function(z) Tab[z,])
        tmp <- lapply(Tab.list,FUN=DRB345.zygosity)
        tmp.df <- lapply(tmp,"[[",1)
        Tab <- as.data.frame(do.call(rbind,tmp.df))
        All.ColNames <- gsub(".1","",colnames(Tab),fixed=T)
        
        #Identify DR345 flagged haplotypes
        DR.Flags <- DRB345.flag(tmp,Tab)
        
        if(Output) {
          if(!is.null(DR.Flags)) {
            colnames(DR.Flags) <- c("Flagged.Locus","Sample.ID") ; rownames(DR.Flags) <- NULL
            Err.Log(Output,"Bad.DRB345.hap") ; cat("\n")
            write.table(DR.Flags,file="Flagged_DRB345_Haplotypes.txt",sep="\t",quote=F,row.names=F,col.names=T)
          }
        }
        cat("\n")
      } else { DRBFLAG <- F }
    } else { DRBFLAG <- NULL }
    
    # Separate locus and allele names if data is formatted as Loci*Allele
    if(HLA==T) { Tab <- CheckLociName(Output,Tab) }
    
  }

# ===================================================================================================================================== ####
# Case-Control Summary ________________________________________________________________________________________________________________ ####
  
  cat(">>>> CASE - CONTROL SUMMARY STATISTICS\n")
  #cat(paste(rep("_",50),collapse=""),"\n")
  if (Trim) { rescall <- paste(Res,"-Field",sep="") } else { rescall <- "Not Defined" }
  Check <- PreCheck(Tab,All.ColNames,rescall,HLA,Verbose,Output)
  if(Output) { write.table(Check,file="Data_Summary.txt",sep=": ",col.names=F,row.names=T,quote=F); rm(Check,rescall) }

# ===================================================================================================================================== ####
# Data Processing and Sanity Checks ___________________________________________________________________________________________________ ####
  
  cat(">>>> DATA PROCESSING AND CHECKS.\n")
  #cat(paste(rep("_",50),collapse=""),"\n")
  
  ## __________________ General processing and checks for any data
  
  # RUN TESTS DEFINITIONS
  if (missing(Run.Tests)) { Run <- c("HWE","H","L","A") } else { Run <- Run.Tests }
  if(!HLA) { 
    if("A" %in% Run) { 
      cat("Not HLA data. Skipping Amino Acid Analysis.\n")
      Run <- Run[-which(Run=="A")]
    }
  }
  
  # BAD Data DEFINITIONS
  if(length(c(which(Tab[,3:ncol(Tab)]==0),which(Tab[,3:ncol(Tab)]==1)))>0) {
    Err.Log(Output,"Bad.Data")
    stop("Analysis Stopped.",call. = F)
  }
  
  # MISSING DATA
  if(Missing == "ignore") {
    cat("Ignoring any missing data.\n")
    cat("Consider setting a missing threshold or running without the haplotype ('H') analysis. A large number of missing data in the haplotype analysis will affect performance, require large amounts of RAM, and in the worst case crash your computer.\n")
    rows.rm <- NULL
  } else {
    if (Missing > 2) {
      if ("H" %in% Run) { cat("The number of allowable missing will affect performance.\nConsider running with a smaller 'Missing' value or without the haplotype ('H') analysis.\ncontinuing......") }
    }
    cat("Removing any missing data. This will affect Hardy-Weinberg Equilibrium test.\n")
    geno.desc <- summaryGeno.2(Tab[,3:ncol(Tab)], miss.val=NAstrings)
    test <- geno.desc[,2] + 2*geno.desc[,3]
    rows.rm <- which(test > Missing)
    if( length(rows.rm) > 0 ) {
      rows.rm <- which(test > Missing)
      ID.rm <- Tab[rows.rm,1]
      Tab <- Tab[-rows.rm,]
      if(Output) { write.table(ID.rm, file="Removed_SampleIDs.txt", sep="\t", row.names=F, col.names=F, quote=F) }
    }
    rm(geno.desc,test,ID.rm)
    if(nrow(Tab)==0) { Err.Log(Output,"TooMany.Missing") ; stop("Analysis Stopped.",call. = F) }
  }
  
  # LOCI SET DEFINITIONS
  if (missing(Loci.Set)) {
    Set <- list(c(3:ncol(Tab)))
  } else {
    Loci.Set <- lapply(Loci.Set,toupper)
    if(CheckLoci(unique(All.ColNames[3:ncol(Tab)]),Loci.Set)$Flag) {
      Err.Log(Output,"Bad.Locus.NA")
      stop("Analysis Stopped.",call. = F)
    } else {
      if(sum(grepl("All",Loci.Set))>0) { Loci.Set[[which(Loci.Set=="All")]] <- unique(All.ColNames[3:ncol(Tab)])  } 
      Set <- lapply(Loci.Set,FUN=function(x)seq(1,ncol(Tab),1)[All.ColNames %in% x]) }
  }
  
  # DATA MERGE AND NUMBER OF LOCI
  if(Output && Merge.Output && All.Pairwise) {
    if(ncol(Tab)>52) {
      cat("You have opted to run all pairwise combinations and merge the final data tables. For a large number of loci, this could take a long time. You have been warned!\n") 
    }
  }
  
  # MULTICORE LIMITATIONS
  if (Cores.Lim!=1L) {
    if(Sys.info()['sysname']=="Windows" & as.numeric(Cores.Lim)>1) { 
      Err.Log(Output,"Cores.Windows")
      stop("Analysis Stopped.",call. = F)
    }
    Cores.Max <- as.integer( floor( parallel::detectCores() * 0.9) )
    if( Cores.Lim > Cores.Max ) { 
      Cores <- Cores.Max
      cat("Adjusting to",Cores.Max,"processor cores.\n")
    } else { Cores <- Cores.Lim }
  } else { Cores <- Cores.Lim }

  ## __________________ HLA specific checks
  if(Trim & !HLA) { cat("Trimming only relevant to HLA data, no trimming performed.\n") }
  if(EVS.rm & !HLA) { cat("Expression variant suffix stripping only relevant to HLA data, no stripping performed.\n") }
  
  if(HLA) {
    
    if(Trim | EVS.rm | "A" %in% Run ) { cat("Running HLA specific functions...\n") }
    
    # Sanity Check for Resoltion if Trim="T" and Trim Data
    if(Trim & CheckHLA(Tab[,3:ncol(Tab)])) {
      cat("--Trimming Data.\n")
      Tab.untrim <- Tab
      Tab[,3:ncol(Tab)] <- apply(Tab[,3:ncol(Tab)],MARGIN=c(1,2),GetField,Res=Res)
      rownames(Tab) <- NULL
    } else if (Trim) {
      Err.Log(Output,"Bad.Format.Trim")
      stop("Analysis Stopped.",call. = F)
    }
    
    # Sanity Check for Expresion Variant Suffix Stripping
    if(EVS.rm & CheckHLA(Tab[,3:ncol(Tab)])) {
      cat("--Stripping Expression Variants Suffixes.\n")
      Tab[,3:ncol(Tab)] <- apply(Tab[,3:ncol(Tab)],MARGIN=c(1,2),gsub,pattern="[[:alpha:]]",replacement="")
      EVS.loci <- as.list(names(EPL))
      EPL <- lapply(EVS.loci,EVSremoval,EPList=EPL)
      names(EPL) <- EVS.loci ; rm(EVS.loci)
    } else if (EVS.rm) {
      Err.Log(Output,"Bad.Format.EVS")
      stop("Analysis Stopped.",call. = F)
    }
    
    if ("A" %in% Run) {
      
      Release <- as.character(unlist(EPL[['Release']]))
      
      # Sanity Check for Known HLA loci
      cat(paste("--Checking loci against ",Release,".\n",sep=""))
      test <- CheckLoci(names(EPL),unique(All.ColNames[3:ncol(Tab)]))
      if( test$Flag ) {
        Err.Log(Output,"Bad.Locus.HLA")
        cat("Problem loci:",test$Loci,"\n")
        stop("Analysis stopped.",call. = F)
      }
      
      # Sanity Check for Known HLA alleles
      cat(paste("--Checking alleles against ",Release,".\n",sep=""))
      test <- CheckAlleles(EPL, Tab[,3:ncol(Tab)], unique(All.ColNames[3:ncol(Tab)]), All.ColNames[3:ncol(Tab)])
      if(sum(unlist(lapply(test,"[[",1)))>0) {
        Err.Log(Output,"Bad.Allele.HLA")
        tmp <- as.character(unlist(lapply(test[which(lapply(test,"[[",1)==T)],"[",2)))
        cat("Problem alleles:",tmp,"\n")
        stop("Analysis stopped.",call. = F)
      }
      
    }
    
  } # End HLA if statement and HLA specific functionalities

  # Multiple Sets And Analysis
  if ( length(Set)>1 & (All.Pairwise | "L" %in% Run | "A" %in% Run ) ) {
    if(!length(Set)==length(unique(unlist(Set)))) { Err.Log(Output,"MultipleSets") }
  }
  

# ===================================================================================================================================== ####
# Write to Parameter File _____________________________________________________________________________________________________________ ####
  
  if(Output) {
    
    if(Data=="HLA_data") { Data.tmp <- "Synthetic Bundled Data Set" } else { Data.tmp <- Data }
    if(HLA && !is.null(DRBFLAG)) { DRB345.tmp <- DRBFLAG } else { DRB345.tmp <- NULL }
    if(HLA) { Trim.tmp <- Trim } else { Trim.tmp <- NULL }
    if(HLA && Trim) { Res.tmp <- Res } else { Res.tmp <- NULL }
    if(HLA) { EVS.rm.tmp <- EVS.rm } else { EVS.rm.tmp <- NULL }
    
    Params.Run <- list(Time = format(Sys.time(), "%a %b %d %X %Y"),
                       BD.Version = as.character(packageVersion("BIGDAWG")),
                       Cores.Used = Cores,
                       File = Data.tmp,
                       Output.Results = Output,
                       Merge = Merge.Output,
                       Return.Object = Return,
                       Display.Results = Verbose,
                       HLA.Data = HLA,
                       DRB345.Parsed = DRB345.tmp,
                       Tests = paste(Run,collapse=","),
                       All.Pairwise = All.Pairwise,
                       Trim = Trim.tmp,
                       Resolution = Res.tmp,
                       Suffix.Stripping = EVS.rm.tmp,
                       Missing.Allowed = Missing,
                       Samples.Removed = length(rows.rm))
    
    Params.Run <- do.call(rbind,Params.Run)
    write.table(Params.Run,file="Run_Parameters.txt",sep=": ", row.names=T, col.names=F, quote=F)
  }

# ===================================================================================================================================== ####
# Hardy Weignberg Equilibrium 'HWE' ___________________________________________________________________________________________________ ####
    
  if ("HWE" %in% Run) {
    
    cat("\n>>>> STARTING HARDY-WEINBERG ANALYSIS...\n")
    #cat(paste(rep("_",50),collapse=""),"\n")
    if(HLA && Trim) { 
      cat("HWE performed at user defined resolution.\n")
    } else { 
      cat("HWE performed at maximum available resolution.\n")
    }
    
    HWE <- HWE.wrapper(Tab,All.ColNames,Output,Verbose)
    BD.out[['HWE']] <- HWE
    rm(HWE)
    
  } #END HARDY-WEINBERG

# ===================================================================================================================================== ####
# Set Loop Begin (loop through each defined locus/loci set) ___________________________________________________________________________ ####
  
  if ( sum( c("H","L","A") %in% Run ) > 0 ) {
  
    cat("\n>>>>>>>>>>>>>>>>>>>>>>>>> Begin Locus Sets <<<<<<<<<<<<<<<<<<<<<<<<<\n\n")
    cat(paste("Your analysis has ", length(Set), " set(s).", sep=""),"\n")
    
    for(k in 1:length(Set)) {
      cat("\n")
      cat(paste(rep(">",35),collapse=""),"Running Set",k,"\n")
      
      cols <- Set[[k]]
      Tabsub <- Tab[,c(1,2,cols)]
      
      #Set Specific Global Variables
      SID <- Tabsub[,1] # sample IDs
      genos <- Tabsub[,3:ncol(Tabsub)] # genotypes
      genos[genos==""] <- NA
      grp <- Tabsub[, 2] # phenotype
      nGrp0 <- length(which(grp==0))*2 #nalleles
      nGrp1 <- length(which(grp==1))*2 #nalleles
      loci <- unique(gsub(".1","",colnames(genos),fixed=T)) # name of loci
      loci.ColNames <- gsub(".1","",colnames(genos),fixed=T) # column names
      nloci <- as.numeric(length(loci)) # number of loci
      SetName <- paste('Set',k,sep="")
      
      if(HLA==T) { genos[genos=='^'] <- "00:00" }
      
      if(Output) {
        
        OutSetDir <- paste(OutDir,"/Set",k,sep="")
        dir.create(OutSetDir)
        setwd(OutSetDir)
        
        Params.set <- list(Set = paste("Set",k),
                       Loci.Run = paste(loci,collapse=","))
        
        Params.set <- do.call(rbind,Params.set)
        write.table(Params.set,file="Set_Parameters.txt",sep=": ", row.names=T, col.names=F, quote=F)
      }
      
      SAFE <- c(ls(),"SAFE")

# ===================================================================================================================================== ####
# Haplotype Analysis 'H' ______________________________________________________________________________________________________________ ####
      
      if ("H" %in% Run) {
        
        #cat(paste(rep("_",50),collapse="","\n"))

        # Sanity check for set length and All.Pairwise=T
        if (nloci<2) {
          Err.Log(Output,"Loci.No")
          stop("Analysis Stopped.", call. = F)
        } else if (All.Pairwise & nloci<=2)  {
          Err.Log(Output,"Loci.No.AP")
          stop("Analysis Stopped.", call. = F) }
        
        Haps.list <- H.MC.wrapper(SID,Tabsub,loci,loci.ColNames,genos,grp,All.Pairwise,Output,Verbose,Cores)
        
        if(All.Pairwise) {
          if(length(BD.out[['H']])>0) { BD.out[['H']] <- c(BD.out[['H']],Haps.list) } else { BD.out[['H']] <- Haps.list }
        } else {
          BD.out[['H']][[SetName]] <- Haps.list
        }
        
        rm(list=ls()[!(ls() %in% SAFE)])
        
      } #END HAPLOTYPE

# ===================================================================================================================================== ####
# Locus Level 'L' _____________________________________________________________________________________________________________________ ####
      
      if ("L" %in% Run) {
        
        #cat(paste(rep("_",50),collapse=""))
        
        L.list <- L.wrapper(nloci,loci,loci.ColNames,genos,grp,nGrp0,nGrp1,Output,Verbose)
        BD.out[['L']][[SetName]] <- list(binned=L.list[['AB']],
                                         freq=L.list[['AF']],
                                         OR=L.list[['OR']],
                                         chisq=L.list[['CS']],
                                         table=L.list[['FB']])
        
        rm(list=ls()[!(ls() %in% SAFE)])
        
      } #END LOCUS

# ===================================================================================================================================== ####
# Amino Acid Level 'A' ________________________________________________________________________________________________________________ ####
      
      if(HLA) {
        if ("A" %in% Run) {
          
          #cat(paste(rep("_",50),collapse=""))

          if(UPL.flag) { cat("Using updated protein exon alignments.\n") }
          
          # Amino Acid Analysis Sanity Checks
          if(Res<2 | !CheckHLA(genos))  {
            cat("You have opted to run the amino acid analysis.\n")
            Err.Log(Output,"Low.Res")
            stop("Analysis stopped.",call. = F)
          }
          
          A.list <- A.wrapper(nloci,loci,loci.ColNames,genos,grp,nGrp0,nGrp1,EPL,Cores,Output,Verbose)
          
          if(Output) {
            ## write to file
            write.table(Release, file = "Set_Parameters.txt", sep="\t", row.names = F, col.names=F, quote = F, append=T)
          }
          
          BD.out[['A']][[SetName]] <- list(log=A.list[['AL']],
                                           binned=A.list[['AB']],
                                           freq=A.list[['AF']],
                                           OR=A.list[['OR']],
                                           chisq=A.list[['CS']],
                                           table=A.list[['FB']])
                      
          rm(list=ls()[!(ls() %in% SAFE)])
          
        } #END AMINO ACID
      }#END if(HLA)

# ===================================================================================================================================== ####
# End Analyses ________________________________________________________________________________________________________________________ ####
      
    }; rm(k)
    
 }# END SET LOOP

  if(Output) {
    
    if(Merge.Output) {
    
      cat("\nMerging data files ...\n")
      if("HWE" %in% Run) { Run <- Run[-which(Run=="HWE")] }
      if( length(Run)>=1 ) { MergeData_Output(BD.out,Run,OutDir) }
      
    }
    
  }
  
# ===================================================================================================================================== ####
  cat("\n>>>>>>>>>>>>>>>>>>>>>>>>>> End Analysis <<<<<<<<<<<<<<<<<<<<<<<<<<\n")
  
  
  if(Output) { setwd(OutDir); save(BD.out, file="Analysis.RData") }
  options(warn=0)
  
  if(Return) { return(BD.out) }
  
}# END FUNCTION


