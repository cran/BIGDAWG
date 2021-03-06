% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/H_support_functions.R
\name{buildHAPsets}
\alias{buildHAPsets}
\title{Haplotype List Builder}
\usage{
buildHAPsets(Combn, genos, loci, loci.ColNames)
}
\arguments{
\item{Combn}{Combination of loci to extraction from genos}

\item{genos}{The genotype columns of the loci set being analyzed.}

\item{loci}{Character vector of unique loci being analyzed.}

\item{loci.ColNames}{Character vector of genos column names.}
}
\description{
Builds table of haplotypes from combinations
}
\note{
This function is for internal BIGDAWG use only.
}
