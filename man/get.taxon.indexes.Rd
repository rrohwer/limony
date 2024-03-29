% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5C_subset_by_taxon.R
\name{get.taxon.indexes}
\alias{get.taxon.indexes}
\title{Get taxon indexes}
\usage{
get.taxon.indexes(my.list, taxa)
}
\arguments{
\item{my.list}{A "flat.list" that includes these names: "names", "av", "sd", "bq", "pd", "br", each of which is an individual matrix.}

\item{taxa}{A character vector of taxa names. Spelling is case sensitive and must be exact.}
}
\value{
A numeric vector of row indexes corresponding to the taxa provided. This can be input into the subset.by.taxa function.
}
\description{
Return a list of row indexes that correspond to a given taxon or vector of taxa.
}
