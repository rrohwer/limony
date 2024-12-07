% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5B_subset_by_sample.R
\name{subset.by.sample}
\alias{subset.by.sample}
\title{Subset by sample (by column index)}
\usage{
\method{subset}{by.sample}(my.list, keep.index)
}
\arguments{
\item{my.list}{Either a "tax.list" or a "flat.list", so a list with these names: these names: "names", "av", "sd", "bq", "pd", "br", each of which contains either a list of taxonomy matrices or a single matrix.}

\item{keep.index}{The index corresponding to sample dates (columns) that you want to keep. Note that keep.index = -index can be used to remove the provided indexes instead.}
}
\value{
Either a "tax.list" or a "flat.list" where all the samples (columns) have been subset to the indexes provided.
}
\description{
Subset a "tax.list" or "flat.list" to the column indexes provided. The function get.sample.indexes can be used to generate the index vector.
}
