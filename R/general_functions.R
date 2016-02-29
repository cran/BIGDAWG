#' Replace absent allele strings
#'
#' Replaces allowable absent allele strings with ^ symbol.
#' @param df Genotypes dataframe.
#' @note This function is for internal BIGDAWG use only.
rmABstrings <- function(df) {
  df[df=="Absent"] <- "^"
  df[df=="absent"] <- "^"
  df[df=="Abs"] <- "^"
  df[df=="ABS"] <- "^"
  df[df=="ab"] <- "^"
  df[df=="Ab"] <- "^"
  df[df=="AB"] <- "^"
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
    binned <- c(NA,NA)
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

