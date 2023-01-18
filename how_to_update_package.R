# RRR
# Run this script to update the package

library(devtools)
library(roxygen2)
library(usethat)
library(usethis)

# Run the first time only to create the package directory
getwd()
setwd("~/Desktop/")
devtools::create("limony")
setwd("limony/")

# re-run every time documentation or code changes
devtools::document()

# re-run to update the included data
limony <- readRDS("~/Desktop/pop/data/limony-IRD/2021-08-25_processing/5A_taxlist_samples_6000.rds")
key <- readRDS("~/Desktop/pop/data/environmental_data/Robin-Refined/seasons/10_limony_package_key.rds")
seasons <- readRDS("~/Desktop/pop/data/environmental_data/Robin-Refined/seasons/10_limony_package_seasons.rds")
usethis::use_data(limony, key, seasons, overwrite = TRUE) # can have multiple objects in this call

# run only the first time to set up the vignette (change doc/*.Rmd to edit and re-install package to test it.)
usethis::use_vignette(name = "introduction")
devtools::build_vignettes(pkg = "limony")

# re-run after pushing changes to github to test the install from github
remove.packages("limony")
devtools::install_github("rrohwer/limony")
library(limony)
data("limony")
# limony:: # to view all script options
devtools::build_vignettes(pkg = "limony")
browseVignettes(package = "limony")
