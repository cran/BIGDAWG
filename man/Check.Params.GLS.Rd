% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_functions.R
\name{Check.Params.GLS}
\alias{Check.Params.GLS}
\title{Check Input Parameters for GLS conversion}
\usage{
Check.Params.GLS(
  Convert,
  File.Output,
  System,
  HZY.Red,
  DRB345.Check,
  Cores.Lim
)
}
\arguments{
\item{Convert}{String Direction for conversion.}

\item{File.Output}{String Type of output.}

\item{System}{String Genetic system (HLA or KIR) of the data being converted}

\item{HZY.Red}{Logical Reduction of homozygote genotypes to single allele.}

\item{DRB345.Check}{Logical Check DR haplotypes for consistency and flag unusual haplotypes.}

\item{Cores.Lim}{Integer How many cores can be used.}
}
\description{
Check input parameters for invalid entries.
}
\note{
This function is for internal use only.
}
