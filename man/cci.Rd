% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat_functions.R
\name{cci}
\alias{cci}
\title{Case-Control Odds ratio calculation and graphing}
\usage{
cci(
  caseexp,
  controlex,
  casenonex,
  controlnonex,
  cctable = NULL,
  graph = TRUE,
  design = "cohort",
  main,
  xlab,
  ylab,
  xaxis,
  yaxis,
  alpha = 0.05,
  fisher.or = FALSE,
  exact.ci.or = TRUE,
  decimal = 2
)
}
\arguments{
\item{caseexp}{Number of cases exposed}

\item{controlex}{Number of controls exposed}

\item{casenonex}{Number of cases not exosed}

\item{controlnonex}{Number of controls not exposed}

\item{cctable}{A 2-by-2 table. If specified, will supercede the outcome and exposure variables}

\item{graph}{If TRUE (default), produces an odds ratio plot}

\item{design}{Specification for graph; can be "case control","case-control", "cohort" or "prospective"}

\item{main}{main title of the graph}

\item{xlab}{label on X axis}

\item{ylab}{label on Y axis}

\item{xaxis}{two categories of exposure in graph}

\item{yaxis}{two categories of outcome in graph}

\item{alpha}{level of significance}

\item{fisher.or}{whether odds ratio should be computed by the exact method}

\item{exact.ci.or}{whether confidence limite of the odds ratio should be computed by the exact method}

\item{decimal}{number of decimal places displayed}
}
\description{
cci function port epicalc version 2.15.1.0 (Virasakdi Chongsuvivatwong, 2012)
}
\note{
This function is for internal BIGDAWG use only.
}
