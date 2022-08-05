# limony
**lake itag measurements over nineteen years**

Please cite: [Rohwer & McMahon. Lake iTag measurements over nineteen years, introducing the limony dataset. _bioRxiv._ 2022](https://www.biorxiv.org/content/10.1101/2022.08.04.502869v1)

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
Or browse the vignette here: [introduction](https://htmlpreview.github.io/?https://github.com/rrohwer/limony/blob/main/introduction.html)

**Contact**:  
Robin.Rohwer@gmail.com  
twitter: @RobinRohwer
