#' Genotype Combination Maker
#'
#' Make data frame of possible genotype combinations
#' @param x Number of alleles.
#' @note This function is for internal BIGDAWG use only.
makeComb <- function(x) {
  tmp <- t(combn(x,2))
  tmp <- rbind(tmp,t(matrix(rep(1:x,2),byrow=T,ncol=x)))
  tmp <- tmp[do.call(order, as.data.frame(tmp)),]
  return(tmp)
}

#' Observed Frequency
#'
#' Get observed frequency of genotypes
#' @param x Single genotype.
#' @param genos.locus Locus genotypes.
#' @note This function is for internal BIGDAWG use only.
getObsFreq <- function(x,genos.locus) {
  
  if(x[1]==x[2]) {
    length(which(genos.locus[,1]==x[1] & genos.locus[,2]==x[2]))
  } else {
    return(sum(length(which(genos.locus[,1]==x[1] & genos.locus[,2]==x[2])),
               length(which(genos.locus[,1]==x[2] & genos.locus[,2]==x[1]))))
  }
}

#' Chi square matrices
#'
#' Chi Square contingency matrix builder with rare cell binning
#' @param Locus Number of alleles.
#' @param genos.sub Genotypes for locus of interest.
#' @param Allele.Freq Allele frequencies.
#' @param Allele.Combn Allele combinations.
#' @note This function is for internal BIGDAWG use only.
getCS.Mat <- function(Locus,genos.sub,Allele.Freq,Allele.Combn) {
  
  nSID <- nrow(genos.sub)
  ColNames <- gsub(".1","",colnames(genos.sub),fixed=T) # column names
  genos.locus <- genos.sub[,which(ColNames==Locus)]
  freq <- Allele.Freq[[Locus]]
  GTYPES <- Allele.Combn[[Locus]]  
  GTYPES <- lapply(seq_len(nrow(GTYPES)), function(x) GTYPES[x,])
  
  #Expected Counts
  freq.Exp <- lapply(GTYPES,FUN=function(x) ifelse(x[1]==x[2],prod(freq[x])*nSID,2*prod(freq[x])*nSID))
  
  #Observed Counts
  GTYPES <- lapply(GTYPES,FUN=function(x) names(freq[x]))
  freq.Obs <- lapply(GTYPES,getObsFreq,genos.locus=genos.locus)
  
  freq.mat <- cbind(do.call(rbind,GTYPES),
                    do.call(rbind,freq.Obs),
                    do.call(rbind,freq.Exp))
  colnames(freq.mat) <- c("Allele.1","Allele.2","Obs","Exp")
  
  #bin rare cells
  freq.bin <- freq.mat[which(as.numeric(freq.mat[,'Exp'])<5),]
  freq.bin <- matrix(data=freq.bin,ncol=ncol(freq.mat),dimnames=dimnames(freq.mat))
  if(nrow(freq.bin)>0) {
    freq.bin <- matrix(c("binned.1","binned.2",sum(as.numeric(freq.bin[,'Obs'])),sum(as.numeric(freq.bin[,'Exp']))),
                       ncol=ncol(freq.bin),
                       dimnames=dimnames(freq.bin))
  }
  
  #Final Matrix for ChiSq
  if(nrow(freq.bin)>0) {
    freq.final <- rbind(freq.mat[which(as.numeric(freq.mat[,'Exp'])>=5),],freq.bin)
  } else {
    freq.final <- freq.mat
  }
  
  #Calculate (Obs - Exp)^2 / Exp
  if(nrow(freq.final)>1) {
    freq.final <- cbind(freq.final,
                        apply(freq.final[,c('Obs','Exp')],
                              MARGIN=1,
                              FUN=function(x) ((as.numeric(x['Obs']) - as.numeric(x['Exp']))^2)/as.numeric(x['Exp'])))
  } else {
    freq.final <- cbind(freq.final,0)
  }
  colnames(freq.final)[ncol(freq.final)] <- "O-E2|E"
  
  return(freq.final)
  
}
