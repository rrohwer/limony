% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5B_subset_by_sample.R
\name{convert.sample.names.to.dates}
\alias{convert.sample.names.to.dates}
\title{Convert sample names to dates}
\usage{
convert.sample.names.to.dates(sample.names)
}
\arguments{
\item{sample.names}{The limony sample names, for ex: colnames(limony$av$kingdom)}
}
\value{
date-formatted sample dates. Note the time zone is GMT-5.
}
\description{
This converts the limony sample names (column names) to a date format.
}
