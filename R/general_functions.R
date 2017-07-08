#' Error Code Display and Logging
#'
#' Displays error codes attributable to data formatting and Locus/Allele naming. Writes to log file.
#' @param Output Logical indicating if Error logging should be written to a file.
#' @param x Log Code.
#' @param y Misc information relevant to error.
#' @note This function is for internal BIGDAWG use only.
Err.Log <- function (Output, x, y=NULL) {
  
  switch(x,
         #Formatting
         Bad.Data = { Error <- "\nYou seem to have data that are 0's or 1's, please considering replacing with another value. Please see vignette." },
         Bad.DRB345.format =  { Error <- "\nYou have included DRB3/4/5 columns, but the alleles calls are not formatted as Locus*Allele. Please see vignette." },
         Bad.DRB345.hap =  { Error <- "\nWe have encountered unanticipated DR haplotypes. Please see the 'Flagged_DRB345_Haplotypes.txt' output file." },
         Uneven.Prefix =  { Error <- "\nIt seems some (not all) of your loci are formatted as Locus*Allele. Please ensure all loci share a similar format." },
         Bad.Format.HLA = { Error <- "\nYour HLA data includes Locus*Allele genotype formatting. Please ensure all known genotypes (including absent calls) follow this format." },
         Bad.Format.Trim = { Error <- "\nYour HLA data does not appear to be formatted properly for trimming. Please see vignette." },
         Bad.Format.EVS = { Error <- "\nYour HLA data does not appear to be formatted properly for EVS stripping. Please see vignette." },
         Case.Con = { Error <- "\nYour data does not appear to contain both cases and controls. Please see vignette." },
         Loci.No = { Error <- "\nYou have opted to run the haplotype analysis with too few loci. Please check Set definitions." },
         Loci.No.AP = { Error <- "\nYou have set All.Pairwise to 'True' but one or more your defined locus sets contain too few loci. Please check Set definitions." },
         Low.Res = { Error <- "\nThe resolution of your HLA data is less than 2 or does not appear to be formatted properly. Please see vignette." },
         High.Res = { Error <- "\nYour HLA does not appear to be formatted properly, >4 fields detected. Please see vignette" },
         #Names
         Bad.Filename = {  Error <- paste("\nBIGDAWG could not locate a file labeled: ",y," in the specificied working directory.",sep="") },
         Bad.Locus.NA = { Error <- "\nYou seem to have specified a locus in the Loci.Set that is not present in your data file." },
         Bad.Locus.HLA = { Error <- "\nThere may be a discrepancy with HLA loci names. Unrecognized locus name(s) encountered." },
         Bad.Allele.HLA = { Error <- "\nThere may be a discrepancy with allele names. Unrecognized allele name(s) encountered." },
         #Other
         MultipleSets = { Error <- "\nWARNING!!! You have opted to run multiple sets with overlapping loci. To avoid duplication of effort and results from the all pairwise haplotype tests, the locus test, and/or the amino acid test(!!!), it is suggested you run these tests separately on either the largest loci set possible or all loci in a given data set." },
         No.Internet = { Error <- "\nYou do not seem to be connected to the internet. CheckRelease() or UpdateRelease() cannot proceed." },
         TooMany.Missing = { Error <- "\nYour data is missing too many values at each locus. Try using Missing='ignore' when running BIGDAWG and avoid haplotype test." }
  )
  cat(Error,"\n")
  if(Output) { write.table(Error,file="Error_Log.txt",sep="\t",quote=F,col.names=F,row.names=F,append=T) }
}

#' Replace absent allele strings
#'
#' Replaces allowable absent allele strings with ^ symbol.
#' @param df Genotypes dataframe.
#' @note This function is for internal BIGDAWG use only.
rmABstrings <- function(df) {
  
  df[,3:ncol(df)] <- apply(df[,3:ncol(df)], MARGIN=c(1,2), FUN=function(x) gsub("ABSENT|Absent|absent|Abs|ABS|ab|Ab|AB|00|00:00|00:00:00|00:00:00:00","^",x) )
  return(df)
  
}

#' Expression Variant Suffix Removal
#'
#' Removes expression variant suffixes from HLA alleles in the exon protein alignment object.
#' @param Locus Locus to be filtered against.
#' @param EPList Exon Protein Alignment Object
#' @note This function is for internal BIGDAWG use only.
EVSremoval <- function(Locus,EPList) {
  if(Locus=='Release') { 
    tmp <- EPList[[Locus]]
    return(tmp)
  } else if(Locus=='RefExons') {
    tmp <- EPList[[Locus]]
    return(tmp)
  } else {
    tmp <- EPList[[Locus]]
    tmp[,'Trimmed'] <- sapply(tmp[,'Trimmed'],gsub,pattern="[[:alpha:]]",replacement="")
    return(tmp)
  }
}

#' DRB345 Column Processing
#'
#' Separates DRB345 column pair into separate columns for each locus
#' @param Tab Data frame of sampleIDs, phenotypes, and genotypes
#' @note This function is for internal BIGDAWG use only.
DRB345.parser <- function(Tab) {
  #Tab Dataset Data-frame

  getCol <- grep("DRB345",colnames(Tab))
  df <- matrix(data="^",nrow=nrow(Tab),ncol=6)
  colnames(df) <- c("DRB3","DRB3.1","DRB4","DRB4.1","DRB5","DRB5.1")
  tmp.1 <- sapply(Tab[,getCol[1]],FUN=GetField,Res=1) ; tmp.2 <- sapply(Tab[,getCol[2]],FUN=GetField,Res=1)
  
  tmp <- list()
  # DRB3
  tmp[[1]] <- unlist(grep("DRB3",Tab[,getCol[1]])) ; tmp[[2]] <- unlist(grep("DRB3",Tab[,getCol[2]]))
  df[tmp[[1]],1] <- Tab[tmp[[1]],getCol[1]] ; df[tmp[[2]],2] <- Tab[tmp[[2]],getCol[2]]
  df[setdiff(1:nrow(df),tmp[[1]]),1] <- "DRB3*^" ; df[setdiff(1:nrow(df),tmp[[2]]),2] <- "DRB3*^"
  df[which(tmp.1=="00"),1] <- paste("DRB3*",Tab[which(tmp.1=="00"),getCol[1]],sep="")
  df[which(tmp.2=="00"),2] <- paste("DRB3*",Tab[which(tmp.2=="00"),getCol[2]],sep="")

  tmp <- list()
  # DRB4
  tmp[[1]] <- unlist(grep("DRB4",Tab[,getCol[1]])) ; tmp[[2]] <- unlist(grep("DRB4",Tab[,getCol[2]]))
  df[tmp[[1]],3] <- Tab[tmp[[1]],getCol[1]] ; df[tmp[[2]],4] <- Tab[tmp[[2]],getCol[2]]
  df[setdiff(1:nrow(df),tmp[[1]]),3] <- "DRB4*^" ; df[setdiff(1:nrow(df),tmp[[2]]),4] <- "DRB4*^"
  df[which(tmp.1=="00"),3] <- paste("DRB4*",Tab[which(tmp.1=="00"),getCol[1]],sep="")
  df[which(tmp.2=="00"),4] <- paste("DRB4*",Tab[which(tmp.2=="00"),getCol[2]],sep="")

  tmp <- list()
  # DRB5
  tmp[[1]] <- unlist(grep("DRB5",Tab[,getCol[1]])) ; tmp[[2]] <- unlist(grep("DRB5",Tab[,getCol[2]]))
  df[tmp[[1]],5] <- Tab[tmp[[1]],getCol[1]] ; df[tmp[[2]],6] <- Tab[tmp[[2]],getCol[2]]
  df[setdiff(1:nrow(df),tmp[[1]]),5] <- "DRB5*^" ; df[setdiff(1:nrow(df),tmp[[2]]),6] <- "DRB5*^"
  df[which(tmp.1=="00"),5] <- paste("DRB5*",Tab[which(tmp.1=="00"),getCol[1]],sep="")
  df[which(tmp.2=="00"),6] <- paste("DRB5*",Tab[which(tmp.2=="00"),getCol[2]],sep="")
  
  # NA's
  df[is.na(Tab[,getCol[1]]),] <- NA ; df[is.na(Tab[,getCol[2]]),] <- NA
  
  Tab.sub <- Tab[,-getCol]
  Tab <- cbind(Tab.sub,df)
  
  return(Tab)
  
}

#' DRB345 haplotype zygosity checker
#'
#' Checks DR haplotypes for correct zygosity and flags unanticipated haplotypes
#' @param x Row of data set data frame following DRB345 parsing
#' @note This function is for internal BIGDAWG use only.
DRB345.zygosity <- function(x) {
  #what about NA or Novels?
  
  Rules <- list("DRB1*01"="^","DRB1*10"="^","DRB1*08"="^",
                "DRB1*03"="DRB3","DRB1*11"="DRB3","DRB1*12"="DRB3","DRB1*13"="DRB3","DRB1*14"="DRB3",
                "DRB1*04"="DRB4","DRB1*07"="DRB4","DRB1*09"="DRB4",
                "DRB1*15"="DRB5","DRB1*16"="DRB5")
  
  
  x <- as.data.frame(x,row.names=rownames(x))
  x.out <- x
  
  x.1F <- apply(x,MARGIN=c(1,2),FUN=GetField,Res=1) # get 1 Field Resolution
  x.1F <- gsub("00","^",x.1F) # substitute "00" with "^"
  
  #DRB1 - get expected DRB3/4/5 genotypes
  DRB1.col <- grep("DRB1",colnames(x.1F))
  
  DRB1.1 <- x.1F[,DRB1.col][1]
  DR.Gtype <- as.character(Rules[DRB1.1])
  
  DRB1.2 <- x.1F[,DRB1.col][2]
  DR.Gtype <- c(DR.Gtype,as.character(Rules[DRB1.2]))
  
  #DRB3 Check
  DRB3.col <- grep("DRB3",colnames(x.1F))
  if( sum(is.na(x.1F[,DRB3.col]))==0 ) {
    
    DRB3.obs <- as.numeric(2 - sum(grepl("\\^",x.1F[,DRB3.col])))
    DRB3.exp <- as.numeric(sum(grepl("DRB3",DR.Gtype)))
    
    A1 <- as.character(x[,DRB3.col[1]]) ; A2 <- as.character(x[,DRB3.col[2]])
    
    if(DRB3.obs!=DRB3.exp) { 
      if(DRB3.obs==2 & DRB3.exp==1 & A1==A2 ) { x.out[,DRB3.col] <- rbind(c("DRB3*^",x[,DRB3.col][1])) ; DR3.flag <- F
      } else { DR3.flag <- T }
    } else { DR3.flag <- F }
    
  } else { DR3.flag <- NA }
  
  #DRB4 Check
  DRB4.col <- grep("DRB4",colnames(x.1F))
  if( sum(is.na(x.1F[,DRB4.col]))==0 ) {
    
    DRB4.obs <- as.numeric(2 - sum(grepl("\\^",x.1F[,DRB4.col])))
    DRB4.exp <- as.numeric(sum(grepl("DRB4",DR.Gtype)))
    
    A1 <- as.character(x[,DRB4.col[1]]) ; A2 <- as.character(x[,DRB4.col[2]])
    
    if(DRB4.obs!=DRB4.exp) { 
      if(DRB4.obs==2 & DRB4.exp==1 & A1==A2 ) { x.out[,DRB4.col] <- rbind(c("DRB4*^",x[,DRB4.col][1])) ; DR4.flag <- F
      } else { DR4.flag <- T }
    } else { DR4.flag <- F }
    
  } else { DR4.flag <- NA }
  
  #DRB5 Check
  DRB5.col <- grep("DRB5",colnames(x.1F))
  if( sum(is.na(x.1F[,DRB5.col]))==0 ) {
    
    DRB5.obs <- as.numeric(2 - sum(grepl("\\^",x.1F[,DRB5.col])))
    DRB5.exp <- as.numeric(sum(grepl("DRB5",DR.Gtype)))
    
    A1 <- as.character(x[,DRB5.col[1]]) ; A2 <- as.character(x[,DRB5.col[2]])
    
    if(DRB5.obs!=DRB5.exp) { 
      if(DRB5.obs==2 & DRB5.exp==1 & A1==A2 ) { x.out[,DRB5.col] <- rbind(c("DRB5*^",x[,DRB5.col][1])) ; DR5.flag <- F
      } else { DR5.flag <- T }
    } else { DR5.flag <- F }
    
  } else { DR5.flag <- NA }
  
  colnames(x.out) <- colnames(x)
  rownames(x.out) <- rownames(x)
  Out.list <- list()
  Out.list[['GTYPE']] <- x.out
  Out.list[['Flag3']] <- DR3.flag
  Out.list[['Flag4']] <- DR4.flag
  Out.list[['Flag5']] <- DR5.flag
  return(Out.list)
  
}

#' DRB345 haplotype zygosity flag check
#'
#' Identify DR345 flagged haplotypes
#' @param tmp Output of DRB345.zygosity
#' @param Tab Data frame of sampleIDs, phenotypes, and genotypes
#' @note This function is for internal BIGDAWG use only.
DRB345.flag <- function(tmp,Tab) {
  
  DR.Flags <- list()
  
  tmp.DR3 <- do.call(rbind,lapply(tmp,"[[",2))
  if(length(which(tmp.DR3==T)>0)) {
    DR.Flags[['DRB3']] <- c("DRB3",paste(Tab[which(tmp.DR3==T),1],collapse=","))
  } else { DR.Flags[['DRB3']] <- NULL }
  
  tmp.DR4 <- do.call(rbind,lapply(tmp,"[[",3))
  if(length(which(tmp.DR4==T)>0)) { 
    DR.Flags[['DRB4']] <- c("DRB4",paste(Tab[which(tmp.DR4==T),1],collapse=",")) 
  } else { DR.Flags[['DRB4']] <- NULL }
  
  tmp.DR5 <- do.call(rbind,lapply(tmp,"[[",4))
  if(length(which(tmp.DR5==T)>0)) { 
    DR.Flags[['DRB5']] <- c("DRB5",paste(Tab[which(tmp.DR5==T),1],collapse=","))
  } else { DR.Flags[['DRB5']] <- NULL }
  
  DR.Flags <- do.call(rbind,DR.Flags)
  return(DR.Flags)
  
}

#' HLA trimming function
#'
#' Trim a properly formatted HLA allele to desired number of fields.
#' @param x HLA allele.
#' @param Res Resolution desired.
#' @note This function is for internal BIGDAWG use only.
GetField <- function(x,Res) {
  Tmp <- unlist(strsplit(as.character(x),":"))
  if (length(Tmp)<2) {
    return(x)
  } else if (Res==1) {
    return(Tmp[1])
  } else if (Res > 1) {
    Out <- paste(Tmp[1:Res],collapse=":")
    return(Out)
  }
}

#' Chi-squared Contingency Table Test
#'
#' Calculates chi-squared contingency table tests and bins rare cells.
#' @param x Contingency table.
#' @note This function is for internal BIGDAWG use only.
RunChiSq <- function(x) {
  
  ### get expected values for cells
  ExpCnts <- chisq.test(as.matrix(x))$expected
  
  ## pull out cells that don't need binning, bin remaining
  #unbinned
  OK.rows <- as.numeric(which(apply(ExpCnts,min,MARGIN=1)>=5))
  if(length(OK.rows)>0) {
    if(length(OK.rows)>=2) {
      unbinned <- x[OK.rows,]
    } else {
      unbinned <- do.call(cbind,as.list(x[OK.rows,]))
      rownames(unbinned) <- rownames(x)[OK.rows]
    }
  } else {
    unbinned <- NULL
  }
  
  #binned
  Rare.rows <- as.numeric(which(apply(ExpCnts,min,MARGIN=1)<5))
  if(length(Rare.rows)>=2) {
    binned <- x[Rare.rows,]
    New.df <- rbind(unbinned,colSums(x[Rare.rows,]))
    rownames(New.df)[nrow(New.df)] <- "binned"
  } else {
    binned <- cbind(NA,NA)
    colnames(binned) <- c("Group.0","Group.1")
    New.df <- x
  }

  if(nrow(New.df)>1) {
  
    # flag if final matrix fails Cochran's rule of thumb (more than 20% of exp cells are less than 5)
    # True = OK ; False = Not good for Chi Square
    ExpCnts <- chisq.test(New.df)$expected
    if(sum(ExpCnts<5)==0){
      flag <- FALSE
    } else if( sum(ExpCnts<5)/sum(ExpCnts>=0)<=0.2 && sum(ExpCnts>=1)>length(ExpCnts) ){
      flag <- FALSE
    } else {
      flag <- TRUE
    }
    
    ## chi square test on binned data
    df.chisq <- chisq.test(New.df)
    Sig <- if(df.chisq$p.value > 0.05) { "NS" } else { "*" }
    
    
    ## show results of overall chi-square analysis
    tmp.chisq <- data.frame(cbind(round(df.chisq$statistic,digits=4),
                                  df.chisq$parameter,
                                  format.pval(df.chisq$p.value),
                                  Sig))
    colnames(tmp.chisq) <- c("X.square", "df", "p.value", "sig")
    
    chisq.out <- list(Matrix = New.df,
                      Binned = binned,
                      Test = tmp.chisq,
                      Flag = flag)
    
    return(chisq.out)
    
  } else {
    
    flag <- TRUE
    tmp.chisq <- data.frame(rbind(rep("NCalc",4)))
    colnames(tmp.chisq) <- c("X.square", "df", "p.value", "sig")
    chisq.out <- list(Matrix = New.df,
                      Binned = binned,
                      Test = tmp.chisq,
                      Flag = flag)
    
  }
  
}

#' Table Maker
#'
#' Table construction of per haplotype for odds ratio, confidence intervals, and pvalues
#' @param x Contingency table with binned rare cells.
#' @note This function is for internal BIGDAWG use only.
TableMaker <- function(x) {
  grp1_sum <- sum(x[,'Group.1'])
  grp0_sum <- sum(x[,'Group.0'])
  grp1_exp <- x[,'Group.1']
  grp0_exp <- x[,'Group.0']
  grp1_nexp <- grp1_sum - grp1_exp
  grp0_nexp <- grp0_sum - grp0_exp
  cclist <- cbind(grp1_exp, grp0_exp, grp1_nexp, grp0_nexp)
  tmp <- as.data.frame(t(cclist))
  names(tmp) <- row.names(x)
  return(tmp)
}

#' Case Control Odds Ratio Calculation from Epicalc
#'
#' Calculates odds ratio and pvalues from 2x2 table
#' @param x List of 2x2 matrices for calculation, output of TableMaker.
#' @note This function is for internal BIGDAWG use only.
cci.pval <- function(x) {
  tmp <- list()
  caseEx <- x[1]
  controlEx <- x[2]
  caseNonEx <- x[3]
  controlNonEx <- x[4]
  table1 <- make2x2(caseEx, controlEx, caseNonEx, controlNonEx)
  tmp1 <- cci(cctable=table1, design = "case-control", graph = FALSE)
  tmp[['OR']] <- round(tmp1$or,digits=2)
  tmp[['CI.L']] <- round(tmp1$ci.or[1],digits=2)
  tmp[['CI.U']] <- round(tmp1$ci.or[2],digits=2)
  tmp[['p.value']] <-  format.pval(chisq.test(table1, correct=F)$p.value)
  tmp[['sig']] <- ifelse(chisq.test(table1, correct=F)$p.value <= 0.05,"*","NS")
  return(tmp)
}

#' Case Control Odds Ratio Calculation from Epicalc list variation
#'
#' Variation of the cci.pvalue function
#' @param x List of 2x2 matrices to apply the cci.pvalue function. List output of TableMaker.
#' @note This function is for internal BIGDAWG use only.
cci.pval.list <- function(x) {
  tmp <- lapply(x, cci.pval)
  tmp <- do.call(rbind,tmp)
  colnames(tmp) <- c("OR","CI.lower","CI.upper","p.value","sig")
  return(tmp)
}

#' Haplotype missing Allele summary function
#'
#' Summary function for identifying missing alleles in a matrix of genotypes.
#' @param geno Matrix of genotypes.
#' @param miss.val Vector of codes for allele missing values.
#' @note This function is for internal BIGDAWG use only and is ported from haplo.stats.
summaryGeno.2 <- function (geno, miss.val = 0) {
  # Ported from R package haplo.stats v 1.7.7
  # Authors: Sinnwell JP, Schaid DJ
  # URL: https://cran.r-project.org/web/packages/haplo.stats/index.html
  n.loci <- ncol(geno)/2
  nr <- nrow(geno)
  geno <- haplo.stats::setupGeno(geno, miss.val)
  loc0 <- numeric(nr)
  loc1 <- numeric(nr)
  loc2 <- numeric(nr)
  for (i in 1:nr) {
    first.indx <- seq(1, (2 * n.loci - 1), by = 2)
    miss.one <- is.na(geno[i, first.indx]) | is.na(geno[i, 
                                                        first.indx + 1])
    miss.two <- is.na(geno[i, first.indx]) & is.na(geno[i, 
                                                        first.indx + 1])
    loc2[i] <- sum(miss.two)
    loc1[i] <- sum(miss.one - miss.two)
    loc0[i] <- sum(!miss.one)
  }
  tbl <- data.frame(missing0 = loc0, missing1 = loc1, missing2 = loc2)
  return(tbl)
}

#' Data Object Merge and Output
#'
#' Whole data set table construction of per haplotype for odds ratio, confidence intervals, and pvalues
#' @param BD.out Output of analysis as list.
#' @param Run Tests that are to be run as defined by Run.Tests.
#' @param OutDir Output directory defined by Results.Dir or default.
#' @note This function is for internal BIGDAWG use only.
MergeData_Output <- function(BD.out,Run,OutDir) {
  
  FM.out <- data.frame(Analysis=character(),
                       Locus=character(),
                       Allele=character(),
                       Group.0=numeric(),
                       Group.1=numeric())
  
  CN.out <- data.frame(Analysis=character(),
                       Locus=character(),
                       Allele=character(),
                       Group.0=numeric(),
                       Group.1=numeric())
  
  OR.out <- data.frame(Analysis=character(),
                       Locus=character(),
                       Allele=character(),
                       OR=numeric(),
                       CI.Lower=numeric(),
                       CI.Upper=numeric(),
                       p.value=numeric(),
                       sig=character())
  
  CS.out <- data.frame(Analysis=character(),
                       Locus=character(),
                       x.square=numeric(),
                       df=numeric(),
                       p.value=numeric(),
                       sig=character())
  
  for(i in Run) {
    
    switch(i,
           H= { TestName <- "Haplotype" },
           L= { TestName <- "Locus" },
           A= { TestName <- "AminoAcid" } )
    
    Test <- BD.out[[i]]
    
    for(k in 1:length(Test)) {
      
      Test.sub <- Test[[k]]
      
      #Frequencies
      tmp <- Test.sub$freq
      if(i=="A") { Allele <- paste(tmp[,'Position'],tmp[,'Residue'],sep="::") }
      switch(i,
             H = { tmp <- cbind(rep(TestName,nrow(tmp)),rep(colnames(tmp)[1],nrow(tmp)),tmp) },
             L = { tmp <- cbind(rep(TestName,nrow(tmp)),tmp) },
             A = { tmp <- cbind(rep(TestName,nrow(tmp)),tmp[,'Locus'],Allele,tmp[,c('Group.0','Group.1')]) })
      colnames(tmp) <- c("Analysis","Locus","Allele","Group.0","Group.1")
      FM.out <- rbind(tmp,FM.out) ; rm(tmp)
      
      #Counts
      tmp <- Test.sub$table
      if(i=="A") { Allele <- paste(tmp[,'Position'],tmp[,'Residue'],sep="::") }
      switch(i,
             H = { tmp <- cbind(rep(TestName,nrow(tmp)),rep(colnames(tmp)[1],nrow(tmp)),tmp) },
             L = { tmp <- cbind(rep(TestName,nrow(tmp)),tmp) },
             A = { tmp <- cbind(rep(TestName,nrow(tmp)),tmp[,'Locus'],Allele,tmp[,c('Group.0','Group.1')]) })
      colnames(tmp) <- c("Analysis","Locus","Allele","Group.0","Group.1")
      CN.out <- rbind(tmp,CN.out) ; rm(tmp)
      
      #Odds Ratios
      tmp <- Test.sub$OR
      if(i=="A") { Allele <- paste(tmp[,'Position'],tmp[,'Residue'],sep="::") }
      switch(i,
             H = { tmp <- cbind(rep(TestName,nrow(tmp)),rep(colnames(tmp)[1],nrow(tmp)),tmp) },
             L = { tmp <- cbind(rep(TestName,nrow(tmp)),tmp) },
             A = { tmp <- cbind(rep(TestName,nrow(tmp)),tmp[,'Locus'],Allele,tmp[,c("OR","CI.lower","CI.upper","p.value","sig")]) })
      colnames(tmp) <- c("Analysis","Locus","Allele","OR","CI.Lower","CI.Upper","p.value","sig")
      OR.out <- rbind(tmp,OR.out) ; rm(tmp)
      
      #ChiSq
      tmp <- Test.sub$chisq
      if(i=="A") { Locus <- paste(tmp[,'Locus'],tmp[,'Position'],sep="::") }
      switch(i,
             H = { tmp <- cbind(rep(TestName,nrow(tmp)), rep(names(Test)[k],nrow(tmp)), tmp) },
             L = { tmp <- cbind(rep(TestName,nrow(tmp)), tmp) },
             A = { tmp <- cbind(rep(TestName,nrow(tmp)), Locus, tmp[,c('X.square','df','p.value','sig')] ) } )
      colnames(tmp)[1:2] <- c("Analysis","Locus")
      rownames(tmp) <- NULL
      CS.out <- rbind(tmp,CS.out); rm(tmp)
      
    }; rm(k)
    
  }; rm(i)
  
  setwd(OutDir)
  
  # Remove redundant entries
  # Especially relevant to multi-set runs
  FM.out <- unique(FM.out)
  CN.out <- unique(CN.out)
  CS.out <- unique(CS.out)
  
  OR.out <- apply(OR.out,MARGIN=c(1,2),as.character)
  OR.out <- unique(OR.out)
  
  write.table(FM.out,file="Merged_Frequencies.txt",sep="\t",col.names=T,row.names=F,quote=F)
  write.table(CN.out,file="Merged_Counts.txt",sep="\t",col.names=T,row.names=F,quote=F)
  write.table(CS.out,file="Merged_ChiSq.txt",sep="\t",col.names=T,row.names=F,quote=F)
  write.table(OR.out,file="Merged_OddsRatio.txt",sep="\t",col.names=T,row.names=F,quote=F)
  
}