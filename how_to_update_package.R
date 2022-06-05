install.packages("devtools")
install.packages("roxygen2")
install.packages("usethat")
library(devtools)
library(roxygen2)
library(usethat)
library(usethis)

getwd()
setwd("~/Desktop/")
devtools::create("limony")
setwd("limony/")
devtools::document() # re-run every time documentation changes

limony <- readRDS("~/Desktop/pop/data/limony-IRD/2021-08-25_processing/5A_taxlist_samples_6000.rds")
usethis::use_data(limony) # can have multiple objects in this call

usethis::use_vignette("introduction")

remove.packages("limony")
devtools::install_github("rrohwer/limony")
library(limony)
data("limony")

