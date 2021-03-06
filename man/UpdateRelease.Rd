% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_wrapper.R
\name{UpdateRelease}
\alias{UpdateRelease}
\title{Update function for protein aligment upon new IMGT HLA data release}
\usage{
UpdateRelease(Force = F, Restore = F, Output = F)
}
\arguments{
\item{Force}{Logical specifiying if update should be forced.}

\item{Restore}{Logical specifying if the original alignment file be restored.}

\item{Output}{Logical indicating if error reporting should be written to file.}
}
\description{
This updates the protein aligment used in checking HLA loci and alleles as well as in the amino acid analysis.
}
