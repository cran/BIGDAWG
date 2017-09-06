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
  
  #Checks for and fixes certain DRB345 errors that are consistent with known DR haplotypes
  
  Rules <- list("DRB1*01"="^","DRB1*10"="^","DRB1*08"="^",
                "DRB1*03"="DRB3","DRB1*11"="DRB3","DRB1*12"="DRB3","DRB1*13"="DRB3","DRB1*14"="DRB3",
                "DRB1*04"="DRB4","DRB1*07"="DRB4","DRB1*09"="DRB4",
                "DRB1*15"="DRB5","DRB1*16"="DRB5")
  
  x.out <- x
  x.1F <- apply(x,MARGIN=c(1,2),FUN=GetField,Res=1) # get 1 Field Resolution
  
  Flag <- NULL
  
  #DRB1 - get expected DRB3/4/5 genotypes
  DRB1.1 <- x.1F[grep("DRB1",colnames(x.1F))[1]]
  DR.Gtype <- as.character(Rules[DRB1.1])
  
  DRB1.2 <- x.1F[grep("DRB1",colnames(x.1F))[2]]
  DR.Gtype <- c(DR.Gtype,as.character(Rules[DRB1.2]))
  
  #DRB3 Check
  DRB3.col <- grep("DRB3",colnames(x.1F))
  DRB3.abs <- grep("\\^",x.1F[DRB3.col])
  if( sum(is.na(x.1F[DRB3.col]))==0 ) {
    
    DRB3.obs <- length(which(sapply(x.1F[DRB3.col],nchar)==7))
    DRB3.exp <- as.numeric(sum(grepl("DRB3",DR.Gtype)))
    A1 <- x.1F[DRB3.col[1]] ; A2 <- x.1F[DRB3.col[2]]
    
    if( DRB3.obs != DRB3.exp ) { 
      if( DRB3.obs==2 && DRB3.exp==1 && A1==A2 ) { x.out[DRB3.col[2]] <- "^"  ; DR3.flag <- F
      } else if( DRB3.obs==2 && DRB3.exp==1 && A1!=A2 ) { DR3.flag <- T
      } else if( DRB3.obs==1 && DRB3.exp==2 ) { x.out[DRB3.col[2]] <- x[DRB3.col[1]]  ; DR3.flag <- F
      } else { DR3.flag <- T }
    } else { DR3.flag <- F }
    
  } else { DR3.flag <- F }
  
  #DRB4 Check
  DRB4.col <- grep("DRB4",colnames(x.1F))
  DRB4.abs <- grep("\\^",x.1F[DRB4.col])
  if( sum(is.na(x.1F[DRB4.col]))==0 ) {
    
    DRB4.obs <- length(which(sapply(x.1F[DRB4.col],nchar)==7))
    DRB4.exp <- as.numeric(sum(grepl("DRB4",DR.Gtype)))
    A1 <- x.1F[DRB4.col[1]] ; A2 <- x.1F[DRB4.col[2]]
    
    if( DRB4.obs != DRB4.exp ) { 
      if( DRB4.obs==2 && DRB4.exp==1 && A1==A2 ) { x.out[DRB4.col[2]] <- "^"  ; DR4.flag <- F
      } else if( DRB4.obs==2 && DRB4.exp==1 && A1!=A2 ) { DR4.flag <- T
      } else if( DRB4.obs==1 && DRB4.exp==2 ) { x.out[DRB4.col[2]] <- x[DRB4.col[1]]  ; DR4.flag <- F
      } else { DR4.flag <- T }
    } else { DR4.flag <- F }
    
  } else { DR4.flag <- F }
  
  #DRB5 Check
  DRB5.col <- grep("DRB5",colnames(x.1F))
  DRB5.abs <- grep("\\^",x.1F[DRB5.col])
  if( sum(is.na(x.1F[DRB5.col]))==0 ) {
    
    DRB5.obs <- length(which(sapply(x.1F[DRB5.col],nchar)==7))
    DRB5.exp <- as.numeric(sum(grepl("DRB5",DR.Gtype)))
    A1 <- x.1F[DRB5.col[1]] ; A2 <- x.1F[DRB5.col[2]]
    
    if( DRB5.obs != DRB5.exp ) { 
      if( DRB5.obs==2 && DRB5.exp==1 && A1==A2 ) { x.out[DRB5.col[2]] <- "^"  ; DR5.flag <- F
      } else if( DRB5.obs==2 && DRB5.exp==1 && A1!=A2 ) { DR5.flag <- T
      } else if( DRB5.obs==1 && DRB5.exp==2 ) { x.out[DRB5.col[2]] <- x[DRB5.col[1]]  ; DR5.flag <- F
      } else { DR5.flag <- T }
    } else { DR5.flag <- F }
  } else { DR5.flag <- F }
  
  # Set Flag
  if(DR3.flag) { Flag <- c(Flag,"DRB3")  }
  if(DR4.flag) { Flag <- c(Flag,"DRB4")  }
  if(DR5.flag) { Flag <- c(Flag,"DRB5")  }
  if(!is.null(Flag)){ names(Flag) <- "DR_Hap_Error" }
  
  # Return Result
  Out <- list()
  colnames(x.out) <- colnames(x)
  rownames(x.out) <- NULL
  Out[['Tab']] <- x.out
  Out[['Flag']] <- ifelse(is.null(Flag),"",paste(Flag,collapse=","))
  return(Out)
  
}

#' DRB345 haplotype zygosity flag check
#'
#' Identify DR345 flagged haplotypes
#' @param tmp Output of DRB345.zygosity
#' @param Tab Data frame of sampleIDs, phenotypes, and genotypes
#' @note This function is for internal BIGDAWG use only.
DRB345.flag <- function(tmp,Tab) {
  
  DR.Flags <- list()
  tmp.Flag <- do.call(rbind,lapply(tmp,"[[",2))
  
  if(sum(sapply(tmp.Flag,grepl,pattern='DRB3'))>0) {
    getDRB.flag <- unlist(sapply(tmp.Flag,grep,pattern='DRB3'))
    DR.Flags[['DRB3']] <- c("DRB3",paste(Tab[getDRB.flag,1],collapse=","))
  } else { DR.Flags[['DRB3']] <- NULL }
  
  if(sum(sapply(tmp.Flag,grepl,pattern='DRB4'))>0) {
    getDRB.flag <- unlist(sapply(tmp.Flag,grep,pattern='DRB4'))
    DR.Flags[['DRB4']] <- c("DRB4",paste(Tab[getDRB.flag,1],collapse=","))
  } else { DR.Flags[['DRB4']] <- NULL }
  
  if(sum(sapply(tmp.Flag,grepl,pattern='DRB5'))>0) {
    getDRB.flag <- unlist(sapply(tmp.Flag,grep,pattern='DRB5'))
    DR.Flags[['DRB5']] <- c("DRB5",paste(Tab[getDRB.flag,1],collapse=","))
  } else { DR.Flags[['DRB5']] <- NULL }
  
  DR.Flags <- do.call(rbind,DR.Flags)
  if(!is.null(DR.Flags)) {
    colnames(DR.Flags) <- c("Flagged.Locus","Sample.ID")
    rownames(DR.Flags) <- NULL
  }
  return(DR.Flags)
  
}
