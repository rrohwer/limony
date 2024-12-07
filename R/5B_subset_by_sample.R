# RRR 2020-11-25
# pull out particular samples
# based on dates or other metadata
# create a list that includes only the subset indices
# or create a list that combines the subset indices
# requires the tax.list and optionally the key (if [DNA] is of interest)

# uses lubridate package

#' Convert sample names to dates
#'
#' This converts the limony sample names (column names) to a date format.
#' @param sample.names The limony sample names, for ex: colnames(limony$av$kingdom)
#' @return date-formatted sample dates. Note the time zone is GMT-5.
#' @export convert.sample.names.to.dates
convert.sample.names.to.dates <- function(sample.names){
  sample.names <- substr(x = sample.names, start = 1, stop = 9)
  sample.names <- lubridate::parse_date_time(x = sample.names, orders = "dmy", tz = "Etc/GMT-5")
  return(sample.names)
}

#' Get sample indexes
#'
#' Get the indexes of samples to either keep or remove. Samples can be selected based on date or based on known method irregularities.
#' @param my.list either a "tax.list" or a "flat.list", so a list with these names: these names: "names", "av", "sd", "bq", "pd", "br", each of which contains either a list of taxonomy matrices or a single matrix.
#' @param start.YY.MM.DD character ymd-formatted date, or "start" to indicate use the earliest date
#' @param end.YY.MM.DD character ymd-formatted date, or "end" to indicate use the latest date
#' @param dates.are.season.range logical. If FALSE, dates used are the range between the start and the end date. If TRUE, dates used are the range each year between the start and end month and day. If TRUE, the year must still be included but can be set to anything.
#' @param remove.prefiltered logical. Remove the non-standard samples that received 10um pre-filtration.
#' @param remove.diff.depths logical. Remove the non-standard samples that were taken at depths different than a 12 m integrated sample.
#' @param remove.diff.loc logical. Remove the non-standard samples that were taken at a different location that the deep hole.
#' @param only.prefiltered logical. Keep only the pre-filtered samples.
#' @param only.diff.depths logical. Keep only the non-standard depths.
#' @param only.diff.loc logical. Keep only the non-standard locations.
#' @param my.key the "key" dataframe supplied with this package.
#' @param remove.low.dna logical. Remove samples that had low-yield DNA extractions.
#' @param low.yield.ng.uL numeric. The yield cutoff for what counts as low (default = 3 ng/uL).
#' @param only.low.dna logical. Keep only low-yield samples.
#' @param selected.seasons character vector. Keep only the listed seasons. Season options are "spring", "clearwater", "early summer", "late summer", "fall", "ice-on".
#' @param selected.invasions character vector. Keep only the listed invasion regimes. Invasion options are "none", "spiny", "zebra".
#' @return A vector of indexes that can be used in the function subset.by.sample()
#' @export get.sample.indexes
get.sample.indexes <- function(my.list,
                           start.YY.MM.DD = "start", end.YY.MM.DD = "end",
                           dates.are.season.range = FALSE,
                           remove.prefiltered = F, remove.diff.depths = F, remove.diff.loc = F,
                           only.prefiltered = F, only.diff.depths = F, only.diff.loc = F,
                           my.key = NA,
                           remove.low.dna = F, low.yield.ng.uL = 3, # samples were diluted to 3ng/uL, not sent if under 1ng/uL
                           only.low.dna = F,
                           selected.seasons = NA, selected.invasions = NA){

  # ---- functions ----

  find.nonstandard.samples <- function(my.samples, remove.prefiltered, remove.diff.depths, remove.diff.loc){
    # grep(pattern = "\\.", x = my.samples, value = T)

    index.pf <- grep(pattern = "pf", x = my.samples, value = F)

    index.D <- grep(pattern = "D", x = my.samples, value = F)
    index.Dec <- grep(pattern = "Dec", x = my.samples, value = F)
    index.DC <- grep(pattern = "DC", x = my.samples, value = F)

    index.D <- setdiff(x = index.D, y = index.Dec) # don't include "Dec" dates, just Depth-labelled
    index.D <- setdiff(x = index.D, y = index.DC) # don't inlude the "DC" samples, which I think is a location?

    index.s <- grep(pattern = "s", x = my.samples, value = F)
    index.s <- union(x = index.s, y = index.DC) # include the "DC" as a spatially different sample

    index.remove <- NULL
    if (remove.prefiltered){ index.remove <- union(index.remove, index.pf) }
    if (remove.diff.depths){ index.remove <- union(index.remove, index.D) }
    if (remove.diff.loc){ index.remove <- union(index.remove, index.s) }

    return(index.remove)
  }

  find.low.dna.samples <- function(my.samples, my.key, dna.ng.ul){

    low.dna <- which(my.key$DNA.ng.uL <= dna.ng.ul)
    low.dna <- my.key$in.R.colnames[low.dna]
    low.dna <- unique(low.dna)

    index.dna <- duplicated(x = c(low.dna,my.samples)) # I think you can do this more elegantly with do.call is.element
    index.dna <- index.dna[-(1:length(low.dna))]
    index.dna <- which(index.dna)

    return(index.dna)
  }

  check.date.format <- function(my.dates, start.YY.MM.DD, end.YY.MM.DD){
    if (length(start.YY.MM.DD) > 1 & !is.numeric(start.YY.MM.DD)){
      return(cat("start date must be character format: \"MM-DD-YY\" \nor the character string \"start\" \nor a numeric vector of pre-selected indexes.\n"))
    }
    if (length(start.YY.MM.DD) > 1){
      return(cat("Using a custom vector of date indexes.\n"))
    }
    if (!is.character(start.YY.MM.DD)){
      return(cat("start date must be character format: \"MM-DD-YY\" \nor the character string \"start\" \nor a numeric vector of pre-selected indexes.\n"))
    }
    if (!is.character(end.YY.MM.DD)){
      return(cat("end date must be character format \"MM-DD-YY\"  \nor set to \"end\" \n"))
    }
    if (!all.equal(order(my.dates), 1:length(my.dates))){
      return(cat("Uh oh... Your dates are not in chronological order, must fix!!\n"))
    }
  }

  find.date.range <- function(my.dates, start.YY.MM.DD, end.YY.MM.DD){
    # my.samples = the in.R.colnames
    # start.YY.MM.DD = either a character date, the word "start", or a numeric vector of date indexes
    # end.YY.MM.DD = either a character date or the word "end" for the last date

    if (length(start.YY.MM.DD) > 1){
      # use if you already found the date indices on your own, but also want to filter by other sample criteria
      return(start.YY.MM.DD)
    }

    if (start.YY.MM.DD == "start"){
      start.YY.MM.DD <- min(my.dates)
      cat("Using earliest date as start date: ", as.character(start.YY.MM.DD), "\n")
    }else{
      start.YY.MM.DD <- lubridate::parse_date_time(x = start.YY.MM.DD, orders = "ymd", tz = "Etc/GMT-5")
    }

    if (end.YY.MM.DD == "end"){
      end.YY.MM.DD <- max(my.dates)
      cat("Using latest date as end date: ", as.character(end.YY.MM.DD), "\n")
    }else{
      end.YY.MM.DD <- lubridate::parse_date_time(x = end.YY.MM.DD, orders = "ymd", tz = "Etc/GMT-5")
    }

    date.range <- lubridate::"%--%"(start.YY.MM.DD, end.YY.MM.DD)
    cat("The date range is: ", as.character(start.YY.MM.DD), "-", as.character(end.YY.MM.DD), "\n")
    index.dates <- which(lubridate::"%within%"(my.dates, date.range))
    return(index.dates)
  }

  find.season.range <- function(my.dates, start.YY.MM.DD, end.YY.MM.DD){
    if (length(start.YY.MM.DD) > 1){
      # use if you already found the date indexes on your own, but also want to filter by other sample criteria
      return(start.YY.MM.DD)
    }

    if (start.YY.MM.DD == "start"){
      start.YY.MM.DD <- paste(min(lubridate::year(my.dates)),1,1,sep = "-")
      start.YY.MM.DD <- lubridate::parse_date_time(x = start.YY.MM.DD, orders = "ymd", tz = "Etc/GMT-5")
    }else{
      start.YY.MM.DD <- lubridate::parse_date_time(x = start.YY.MM.DD, orders = "ymd", tz = "Etc/GMT-5")
    }

    if (end.YY.MM.DD == "end"){
      end.YY.MM.DD <- paste(max(lubridate::year(my.dates)),12,31,sep = "-")
      end.YY.MM.DD <- lubridate::parse_date_time(x = end.YY.MM.DD, orders = "ymd", tz = "Etc/GMT-5")
    }else{
      end.YY.MM.DD <- lubridate::parse_date_time(x = end.YY.MM.DD, orders = "ymd", tz = "Etc/GMT-5")
    }

    my.yrs <- unique(lubridate::year(my.dates))
    my.start.month <- lubridate::month(start.YY.MM.DD)
    my.start.day <- lubridate::day(start.YY.MM.DD)
    my.end.month <- lubridate::month(end.YY.MM.DD)
    my.end.day <- lubridate::day(end.YY.MM.DD)

    index.dates <- NULL
    for(y in my.yrs){
      season.start <- paste(y, my.start.month, my.start.day, sep = "-")
      season.end <- paste(y, my.end.month, my.end.day, sep = "-")
      season.index <- find.date.range(my.dates = my.dates, start.YY.MM.DD = season.start, end.YY.MM.DD = season.end)
      index.dates <- c(index.dates, season.index)
    }
    return(index.dates)
  }

  find.phenology.seasons <- function(my.key, selected.seasons){
    index <- my.key$Season %in% selected.seasons
    return(which(index))
  }

  find.invasion.periods <- function(my.key, selected.invasions){
    index <- my.key$Invasion %in% selected.invasions
    return(which(index))
  }

  # ---- actions ----

  # my.list can be tax.list or flat.list
  if ( !is.matrix(my.list$bq) ){
    my.samples <- colnames(my.list$bq$Kingdom)
  }else{
    my.samples <- colnames(my.list$bq)
  }

  my.dates <- convert.sample.names.to.dates(sample.names = my.samples)

  if (remove.prefiltered | remove.diff.depths | remove.diff.loc){
    index.nonstandard.toss <- find.nonstandard.samples(my.samples = my.samples, remove.prefiltered = remove.prefiltered, remove.diff.depths = remove.diff.depths, remove.diff.loc = remove.diff.loc)
  }else{
    index.nonstandard.toss <- NULL
  }

  if (remove.low.dna){
    if (is.na(my.key)[1]){
      return(cat("Need to supply the sample key table (my.key = ) in order to filter by extraction yield.\n"))
    }
    index.lowyield.toss <- find.low.dna.samples(my.samples = my.samples, my.key = my.key, dna.ng.ul = low.yield.ng.uL)
  }else{
    index.lowyield.toss <- NULL
  }

  if (!is.na(selected.seasons[1])){
    if (is.na(my.key)[1]){
      return(cat("Need to supply the sample key table (my.key = ) in order to filter by season.\n"))
    }
    season.options <- unique(my.key$Season)
    if(sum(duplicated(c(season.options, selected.seasons))) < length(selected.seasons)){
      return(cat("You spelled a season wrong. Correct seasons are:\n",paste(season.options,collapse = ", "),"\n", "You supplied:\n",paste(selected.seasons, collapse = ", "),"\n"))
    }
    index.seasons <- find.phenology.seasons(my.key = my.key, selected.seasons = selected.seasons)
  }else{
    index.seasons <- NA
  }

  if (!is.na(selected.invasions[1])){
    if (is.na(my.key)[1]){
      return(cat("Need to supply the sample key table (my.key = ) in order to filter by species invasion\n"))
    }
    invasion.options <- unique(my.key$Invasion)
    if(sum(duplicated(c(invasion.options, selected.invasions))) < length(selected.invasions)){
      return(cat("You spelled an invasion regime wrong. Correct invasions are:\n",paste(invasion.options,collapse = ", "),"\n", "You supplied:\n",paste(selected.invasions, collapse = ", "),"\n"))
    }
    index.invasions <- find.invasion.periods(my.key = my.key, selected.invasions = selected.invasions)
  }else{
    index.invasions <- NA
  }

  if (!dates.are.season.range){
    index.daterange <- find.date.range(my.dates = my.dates, start.YY.MM.DD = start.YY.MM.DD, end.YY.MM.DD = end.YY.MM.DD)
  }else{
    index.daterange <- find.season.range(my.dates = my.dates, start.YY.MM.DD = start.YY.MM.DD, end.YY.MM.DD = end.YY.MM.DD)
  }

  index.keep <- lubridate::setdiff(x = index.daterange, y = index.nonstandard.toss)
  index.keep <- lubridate::setdiff(x = index.keep, y = index.lowyield.toss)

  if (only.prefiltered | only.diff.depths | only.diff.loc){
    index.nonstandard.keep <- find.nonstandard.samples(my.samples = my.samples, remove.prefiltered = only.prefiltered, remove.diff.depths = only.diff.depths, remove.diff.loc = only.diff.loc)
    index.keep <- lubridate::intersect(x = index.keep, y = index.nonstandard.keep)
  }

  if (only.low.dna){
    index.lowyield.keep <- find.low.dna.samples(my.samples = my.samples, my.key = my.key, dna.ng.ul = low.yield.ng.uL)
    index.keep <- lubridate::intersect(x = index.keep, y = index.lowyield.keep)
  }

  if (!is.na(index.invasions[1])){
    index.keep <- lubridate::intersect(x = index.keep, y = index.invasions)
  }
  if (!is.na(index.seasons[1])){
    index.keep <- lubridate::intersect(x = index.keep, y = index.seasons)
  }

  cat("Returning the column indexes for", length(index.keep), "total samples that you want to keep.\n")

  return(index.keep)
}

#' Subset by sample (by column index)
#'
#' Subset a "tax.list" or "flat.list" to the column indexes provided. The function get.sample.indexes can be used to generate the index vector.
#' @param my.list Either a "tax.list" or a "flat.list", so a list with these names: these names: "names", "av", "sd", "bq", "pd", "br", each of which contains either a list of taxonomy matrices or a single matrix.
#' @param keep.index The index corresponding to sample dates (columns) that you want to keep. Note that keep.index = -index can be used to remove the provided indexes instead.
#' @return Either a "tax.list" or a "flat.list" where all the samples (columns) have been subset to the indexes provided.
#' @export subset.by.sample
subset.by.sample <- function(my.list, keep.index){
  # my.list can be tax.list or flat.list or y.vals (flat.list + extra y.max element)
  if ( !is.matrix(my.list$bq) ){
    names(my.list)
    for (t in 1:length(my.list$bq)){
      my.list$av[[t]] <- my.list$av[[t]][ ,keep.index, drop = F]
      my.list$sd[[t]] <- my.list$sd[[t]][ ,keep.index, drop = F]
      my.list$bq[[t]] <- my.list$bq[[t]][ ,keep.index, drop = F]
      my.list$pd[[t]] <- my.list$pd[[t]][ ,keep.index, drop = F]
      my.list$br[[t]] <- my.list$br[[t]][ ,keep.index, drop = F]
    }
  }else{
    my.list$av <- my.list$av[ ,keep.index, drop = F]
    my.list$sd <- my.list$sd[ ,keep.index, drop = F]
    my.list$bq <- my.list$bq[ ,keep.index, drop = F]
    my.list$pd <- my.list$pd[ ,keep.index, drop = F]
    my.list$br <- my.list$br[ ,keep.index, drop = F]
    if (any(duplicated(c(names(my.list), "filled.sd")))){ # this element exists in y.vals data structure
      my.list$filled.sd <- my.list$filled.sd[ ,keep.index, drop = F]
    }
  }
  return(my.list)
}

#' Group by sample
#'
#' Lump samples together into a single group. For example, you could group all June samples into one sample.
#'
#' @param my.list Either a "tax.list" or a "flat.list."
#' @param group.index.list A list where each element is a vector of sample (column) indexes that should be grouped together.
#' @param new.colname.vect A vector of new column names. The length of this vector should match the length of the group.index.list.
#' @return Either a "tax.list" or a "flat.list" where the the samples have been grouped together. Taxon abundances are summed, standard deviation is propagated, and quantification limits are recalculated.
#' @export group.by.sample
group.by.sample <- function(my.list, group.index.list, new.colname.vect){

  propagate.error.sums <- function(error.vector){
    # note: there's a package called propagate that can figure out more complex formulas. seems easier to do by hand since formulas simple.
    sx <- sqrt( sum( error.vector ^ 2 , na.rm = T) )
    return(sx)
  }

  num.groups <- length(group.index.list)
  new.list <- my.list # new list is same format, but same number columns as groups (will re-fill and re-name the columns)

  if ( !is.matrix(my.list$bq) ){
    names(my.list)
    for (t in 1:length(my.list$bq)){
      new.list$av[[t]] <- new.list$av[[t]][ ,1:num.groups, drop = F]
      new.list$sd[[t]] <- new.list$sd[[t]][ ,1:num.groups, drop = F]
      new.list$bq[[t]] <- new.list$bq[[t]][ ,1:num.groups, drop = F]
      new.list$pd[[t]] <- new.list$pd[[t]][ ,1:num.groups, drop = F]
      new.list$br[[t]] <- new.list$br[[t]][ ,1:num.groups, drop = F]
      colnames(new.list$av[[t]]) <- new.colname.vect
      colnames(new.list$sd[[t]]) <- new.colname.vect
      colnames(new.list$bq[[t]]) <- new.colname.vect
      colnames(new.list$pd[[t]]) <- new.colname.vect
      colnames(new.list$br[[t]]) <- new.colname.vect
    }
    for(g in 1:num.groups){
      temp.list <- subset.by.sample(my.list = my.list, keep.index = group.index.list[[g]])
      for (t in 1:length(temp.list$bq)){
        new.list$av[[t]][ ,g] <- apply(X = temp.list$av[[t]], MARGIN = 1, FUN = mean)

        sd.of.avs <- apply(X = temp.list$av[[t]], MARGIN = 1, FUN = sd)
        sd.propagated <- apply(X = temp.list$sd[[t]], MARGIN = 1, FUN = propagate.error.sums)
        sd.propagated <- sd.propagated / nrow(temp.list$sd[[t]]) # mean is sums divided by a constant n
        index <- sd.propagated > sd.of.avs
        cat("prop is bigger",sum(index)/nrow(temp.list$sd[[t]]),"\n")
        larger.sd <- sd.of.avs
        larger.sd[index] <- sd.propagated[index]
        new.list$sd[[t]][ ,g] <- larger.sd

        new.list$bq[[t]][ ,g] <- apply(X = temp.list$bq[[t]], MARGIN = 1, FUN = all)
        new.list$pd[[t]][ ,g] <- apply(X = temp.list$pd[[t]], MARGIN = 1, FUN = all)
        new.list$br[[t]][ ,g] <- apply(X = temp.list$br[[t]], MARGIN = 1, FUN = all)
      }
    }
  }else{
    new.list$av <- new.list$av[ ,1:num.groups, drop = F]
    new.list$sd <- new.list$sd[ ,1:num.groups, drop = F]
    new.list$bq <- new.list$bq[ ,1:num.groups, drop = F]
    new.list$pd <- new.list$pd[ ,1:num.groups, drop = F]
    new.list$br <- new.list$br[ ,1:num.groups, drop = F]
    colnames(new.list$av) <- new.colname.vect
    colnames(new.list$sd) <- new.colname.vect
    colnames(new.list$bq) <- new.colname.vect
    colnames(new.list$pd) <- new.colname.vect
    colnames(new.list$br) <- new.colname.vect
    if (any(duplicated(c(names(my.list), "filled.sd")))){ # this element exists in y.vals data structure
      new.list$filled.sd <- new.list$filled.sd[ ,1:num.groups, drop = F]
      colnames(new.list$filled.sd) <- new.colname.vect
    }
    for(g in 1:num.groups){
      temp.list <- subset.by.sample(my.list = my.list, keep.index = group.index.list[[g]])

      new.list$av[ ,g] <- apply(X = temp.list$av, MARGIN = 1, FUN = mean)

      sd.of.avs <- apply(X = temp.list$av, MARGIN = 1, FUN = sd)
      sd.propagated <- apply(X = temp.list$sd, MARGIN = 1, FUN = propagate.error.sums)
      sd.propagated <- sd.propagated / nrow(temp.list$sd) # mean is sums divided by a constant n (prop sums / const)
      index <- sd.propagated > sd.of.avs
      # cat("prop is bigger",sum(index)/nrow(temp.list$sd),"\n") # propagated error is only bigger at kingdom level, it seems
      larger.sd <- sd.of.avs
      larger.sd[index] <- sd.propagated[index]
      new.list$sd[ ,g] <- larger.sd

      new.list$bq[ ,g] <- apply(X = temp.list$bq, MARGIN = 1, FUN = all)
      new.list$pd[ ,g] <- apply(X = temp.list$pd, MARGIN = 1, FUN = all)
      new.list$br[ ,g] <- apply(X = temp.list$br, MARGIN = 1, FUN = all)

      if (any(duplicated(c(names(my.list), "filled.sd")))){ # this element exists in y.vals data structure
        new.list$filled.sd[ ,g] <- new.list$sd[ ,g] # there are no NA's in the SDs anymore b/c takes the sd.of.avs if a prop sd doesn't exist...I think,didn't test yet
      }
    }
  }
  return(new.list)
}



