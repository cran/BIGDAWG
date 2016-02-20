#' Haplotype Analysis Function
#'
#' This is the workhorse function for the haplotype analysis.
#' @param genos.sub The genotype columns of the loci(locus) set being analyzed. 
#' @param grp Case/Control or Phenotype groupings.
#' @note This function is for internal BIGDAWG use only.
H <- function(genos.sub,grp) {
  
    cat("Estimating Haplotypes ...\n")
    
    loci.sub <- unique(gsub(".1","",colnames(genos.sub),fixed=T))
    nloci.sub <- as.numeric(length(loci.sub))
    Haplotype <- paste(loci.sub,collapse="~")  
    
    ### estimate haplotypes
    Tab.out <- haplo.stats::haplo.group(group = grp, genos.sub, locus.label = loci.sub)
    
    ## extract haplotype freqs for cases and controls
    rescol <- c((nloci.sub + 2), (nloci.sub + 3)) 
    haps <- (Tab.out$group.df[,rescol])
    haps[is.na(haps)] <- 0
    names(haps) <- c("Group.0", "Group.1")
    hapres <- Tab.out$group.df
    loci_haps <- c(1:nloci.sub)
    row.names(haps) <- do.call("paste", c(hapres[loci_haps], sep = "~"))
    hapsrnd <- round(haps, 5)
    
    ## convert freqs to counts
    go_cnt <- (data.frame(Tab.out$group.count))[1, 2]
    g1_cnt <- (data.frame(Tab.out$group.count))[2, 2]
    haps_counts <- haps
    haps_counts[,'Group.0'] <- round(2 * go_cnt * (haps_counts[,'Group.0']))
    haps_counts[,'Group.1'] <- round(2 * g1_cnt * (haps_counts[,'Group.1']))
    
    ### get expected values for cells, bin small cells, and run chi square
    Result <- RunChiSq(haps_counts)
    haps_binned <- Result$Binned
    Final_binned <- Result$Matrix
    overall.chisq <- Result$Test
    
    ## make a nice table of ORs, ci, p values
    ccdat <-TableMaker(Final_binned)
    ORout <- lapply(ccdat, cci.pval) #OR
    ORout <- do.call(rbind,ORout)
    
    ## fix tables so that haplotype row names become a separate column
    #haps_binned
    haps_binned_fix <- cbind(rownames(haps_binned),haps_binned)
    colnames(haps_binned_fix) <- c(Haplotype,colnames(haps_binned))
    rownames(haps_binned_fix) <- NULL
    
    #final_binned
    Final_binned_fix <- cbind(rownames(Final_binned),Final_binned)
    colnames(Final_binned_fix) <- c(Haplotype,colnames(Final_binned))
    rownames(Final_binned_fix) <- NULL
    
    #hapsrnd
    hapsrnd_fix <- cbind(rownames(hapsrnd),hapsrnd)
    colnames(hapsrnd_fix) <- c(Haplotype,colnames(hapsrnd))
    rownames(hapsrnd_fix) <- NULL
    
    #ORout
    ORout_fix <- cbind(rownames(ORout),ORout)
    colnames(ORout_fix) <- c(Haplotype,colnames(ORout))
    rownames(ORout_fix) <- NULL
    
    H.tmp <- list()
    H.tmp[['freq']] <- hapsrnd_fix # rounded frequencies  
    H.tmp[['binned']] <- haps_binned_fix # binned haplotypes
    H.tmp[['OR']] <- ORout_fix # odd ratio table
    H.tmp[['chisq']] <- overall.chisq # chi sq test statistic
    H.tmp[['table']] <- Final_binned_fix # final table for chi sq

    return(H.tmp)
  
}
