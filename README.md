# limony
lake itag measurements over nineteen years

Please cite: INSERT WHEN HAVE

This repo contains scripts that were used to process the limony dataset, and also a package. Loading the package provides both the limony data and functions for subsetting its nested list structure. 

To load the package, first make sure you have the pre-requisites installed:
```
install.packages("devtools")
install.packages("lubridate")
```

Now install the limony package like so:
```
devtools::install_github("rrohwer/limony", build_vignettes = TRUE)
```

Find the limony dataset and directions with:
```
library(limony)
data("limony")
browseVignettes(package = "limony")
```
Or find the vignette here: ugh where is it
