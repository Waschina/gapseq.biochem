% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Reaction.R
\name{getRxnEquation}
\alias{getRxnEquation}
\title{Get reaction equations}
\usage{
getRxnEquation(
  db,
  rxn,
  highlight = NULL,
  use.ids = FALSE,
  format.style = FALSE
)
}
\arguments{
\item{db}{Database of class \link{gapseqDB-class}}

\item{rxn}{Reaction ID(s)}

\item{highlight}{Compound IDs whose names/IDs should be highlighted in red.
(only works if 'format.style' is TRUE)}

\item{use.ids}{Use compound IDs instead of compound names?}

\item{format.style}{Should the output be formatted (highlights and grey
background for external metabolites)?}
}
\description{
Nicely formatted reaction equations with compound names or IDs.
}
