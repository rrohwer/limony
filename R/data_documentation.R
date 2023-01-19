# RRR
# documentation for the data objects

#' limony
#'
#' The nested list structure that contains all OTU tables.
#' For a visual diagram, see Figure 3 in the limony data announcement.
#' Please cite the limony data announcement if you use this data:
#' Rohwer & McMahon. Lake iTag measurements over nineteen years, introducing the limony dataset. bioRxiv. 2022
#' @format ## `limony`
#' A list with 6 elements:
#' \describe{
#'    \item{names}{Taxonomic Assignments (character matrices)}
#'    \item{av}{average abundance (integer matrices)}
#'    \item{sd}{standard deviation (integer matrices)}
#'    \item{bq}{below quantification (boolean matrices)}
#'    \item{pd}{partial detection (boolean matrices)}
#'    \item{br}{below resolution (boolean matrices)}
#' }
#' @source Rohwer & McMahon. Lake iTag measurements over nineteen years, introducing the limony dataset. bioRxiv. 2022
"limony"

#' key
#'
#' A data frame that includes information about each sample
#' If you use the season information, please cite the paper that defined the seasons:
#' Rohwer et al. Species invasions shift microbial phenology in a two-decade freshwater time series. bioRxiv. 2022
#' @format ## `key`
#' \describe{
#'    \item{in.R.colnames}{exact matches to the sample names in limony matrix columns}
#'    \item{in.R.order}{the order of columns in the limony matrices}
#'    \item{Year}{Sample year}
#'    \item{Month}{Sample month}
#'    \item{Day}{Sample day}
#'    \item{Date}{Sample date, with time zone summer in Wisconsin}
#'    \item{Invasion}{Invasion regime}
#'    \item{Season}{Season as defined in Rohwer et al. 2023}
#'    \item{DNA.ng.uL}{Concentration of DNA extraction (can be used to filter low-yield samples)}
#'    \item{Rep.Type}{Type of replicate}
#'    \item{Color.Invasion}{Colors used in Rohwer et al. 2023}
#'    \item{Color.Season}{Colors used in Rohwer et al. 2023}
#'    \item{yDay}{Julian day of year}
#'    \item{Days.Since.Spring}{Days since the start of Spring}
#'    \item{Days.Since.Clearwater}{Days since the start of Clearwater}
#'    \item{Days.Since.Early.Summer}{Days since the start of Early Summer}
#'    \item{Days.Since.Late.Summer}{Days since the start of Late Summer}
#'    \item{Days.Since.Fall}{Days since the start of Fall}
#'    \item{Days.Since.Ice.On}{Days since the start of the previous Ice-On}
#'}
#' @source Rohwer et al. Species invasions shift microbial phenology in a two-decade freshwater time series. bioRxiv. 2022
"key"

#' seasons
#'
#' A data frame that includes start dates for each season each year
#' If you use the season information, please cite the paper that defined the seasons:
#' Rohwer et al. Species invasions shift microbial phenology in a two-decade freshwater time series. bioRxiv. 2022
#' @format ## `seasons`
#' \describe{
#'    \item{Year}{The observation year}
#'    \item{Spring}{Start dates for the season Spring}
#'    \item{Clearwater}{Start dates for the season Clearwater}
#'    \item{Early.Summer}{Start dates for the season Early Summer}
#'    \item{Late.Summer}{Start dates for the season Late Summer}
#'    \item{Fall}{Start dates for the season Fall}
#'    \item{Ice.On}{Start dates for the season Ice-on}
#'}
#' @source Rohwer et al. Species invasions shift microbial phenology in a two-decade freshwater time series. bioRxiv. 2022
"seasons"
