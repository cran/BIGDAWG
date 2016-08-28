#' Haplotype Wrapper
#'
#' Wrapper for main H function
#' @param SID Character vector of subject IDs.
#' @param Tabsub Data frame of genotype calls for set being analyzed.
#' @param loci Character vector of unique loci being analyzed.
#' @param loci.ColNames Character vector of genos column names.
#' @param genos The genotype columns of the loci set being analyzed. 
#' @param grp Case/Control or Phenotype groupings.
#' @param All.Pairwise Haplotype argument carryover from main BIGDAWG function
#' @param Output Data return carryover from main BIGDAWG function
#' @param Verbose Summary display carryover from main BIGDAWG function
#' @note This function is for internal BIGDAWG use only.
H.wrapper <- function(SID,Tabsub,loci,loci.ColNames,genos,grp,All.Pairwise,Output,Verbose) {
  
  HAPsets <- list() ; Haps.list <- list()
  HAPsets[[paste(loci,collapse='~')]] <- genos
  
  # Define Pairwise Combinations to Run When Selected
  if(All.Pairwise) {
    
    # Define Combinations
    Combos <- t(combn(loci,2))
    
    cat("\nYou have opted to include all pairwise combinations for the haplotype analysis.\n")
    cat("There are", nrow(Combos), "possible locus combinations to run.\n" )
    
    # Define Pairwise Sets
    for(s in 1:nrow(Combos)) {
      Set.H <- loci.ColNames %in% Combos[s,]
      HAPsets[[paste(Combos[s,1],Combos[s,2],sep="~")]] <- genos[,c(Set.H)]
    }; rm(s)
    
  }
  
  for(h in 1:length(HAPsets)) {
    
    # Get haplotype loci set for analysis
    genos.sub <- HAPsets[[h]]
    
    # Run Analysis
    H.list <- H(genos.sub,grp)
    
    if(Output) {
      
      # File names for output
      name1 <- paste("haplotype_freqs.",names(HAPsets)[h],".txt",sep="")
      name2 <- paste("haplotype_binned.",names(HAPsets)[h],".txt",sep="")
      name3 <- paste("haplotype_OR.",names(HAPsets)[h],".txt",sep="")
      name4 <- paste("haplotype_chisq.",names(HAPsets)[h],".txt",sep="")    
      name5 <- paste("haplotype_table.",names(HAPsets)[h],".txt",sep="") 
      name6 <- paste("haplotype_bySubject.",names(HAPsets)[h],".txt",sep="") 
      
      ## write to file
      write.table(H.list[['freq']], name1, sep="\t", quote = F, row.names=F, col.names=T)
      write.table(H.list[['binned']], name2, sep="\t", quote = F, row.names=F, col.names=T)
      write.table(H.list[['OR']], name3, sep="\t", quote = F, row.names=F, col.names=T)
      write.table(H.list[['chisq']], name4, sep="\t", row.names = F, quote = F)
      write.table(H.list[['table']], name5, sep="\t", row.names = F, quote = F)
      
      Haplotypes <- cbind(SID,H.list[['Haplotypes']])
      colnames(Haplotypes)[1] <- colnames(Tabsub)[1]
      write.table(Haplotypes, name6, sep="\t", row.names = F, quote = F)
      
    }
    
    cat("> HAPLOTYPE ANALYSIS COMPLETED:",names(HAPsets)[h],"\n")
    if(Verbose) {
      overall.chisq <- H.list[['chisq']]
      overall.chisq$X.square <- round(as.numeric(levels(overall.chisq$X.square)),digits=5)
      print(overall.chisq, row.names=F)
      cat("\n")
    }
    
    SetName <- gsub(" ","",names(HAPsets)[h])
    Haps.list[[SetName]] <- H.list
    
  }# END hapset Loop
  
  return(Haps.list)
  
}