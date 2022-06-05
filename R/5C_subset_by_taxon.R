# RRR 5/18/20
# Pull out particular taxa (can provide multiple names)
# Renormalize or not
# also provides functions to find and walk through taxon names (lists must be flattened first)

# Use functions in this order
# 1. flatten taxa lists to only have seqID level
# 2. get indexes to keep
# 3. subset to those indexes
# 4. re-group (with the separate grouping script)

#' Flatten a "tax.list" to a "flat.list"
#' 
#' Take a nested list with every taxonomy level and flatten it into a single list with only 1 taxonomy level.
#' @param my.list A "tax.list" that has these names: "names", "av", "sd", "bq", "pd", "br", each of which is a nested list that has a matrix for each taxonomy level.
#' @param finest.level Either a numeric taxonomy level or a character taxon name. This defines the taxonomy level that the list is flattened to. By default, the OTU-level is chosen.
#' @return A "flat.list" that includes these names: "names", "av", "sd", "bq", "pd", "br", each of which is an individual matrix at the chosen taxonomy level.
#' @export
flatten.to.single.level <- function(my.list, finest.level = length(my.list$bq)){
  if (is.character(finest.level)){
    finest.level <- get.taxon.level(my.list = my.list, taxon = finest.level)
  }
  
  for (e in 1:length(my.list)){
    my.list[[e]] <- my.list[[e]][[finest.level]]
  }
  return(my.list)
}

#' Get the taxonomy level of a given taxon name
#' 
#' Get the numeric level of the provided taxon name. Kingdom = 1, Phylum = 2, Class = 3, Order = 4, Family/Lineage = 5, Genus/Clade = 6, Species/Tribe = 7, SeqID = 8.
#' @param my.list Either a "tax.list" or a "flat.list"
#' @param taxon An exactly spelled taxon name. This is case sensitive.
#' @return The numeric taxon level.
#' @export
get.taxon.level <- function(my.list, taxon){
  # my.list can be tax.list or flat.list
  # this will not search through OTU names
  
  if (!is.matrix(my.list$names)){
    my.names <- my.list$names$`Species/Tribe`
  }else{
    my.names <- my.list$names
    if (ncol(my.names) > 7){
      my.names <- my.names[ ,1:7]
    }
  }
  
  index <- NA
  
  for (t in 1:(ncol(my.names))){
    index <- which(my.names[ ,t] == taxon)
    if (length(index) > 0){
      taxa.level <- t
      break
    }
  }
  if (length(index) < 1){
    return(cat(taxon, "is not a taxon name. These names include the beginning of what you typed:\n\n", 
               grep(pattern = substr(x = taxon, start = 1, stop = 5), x = unique(as.vector(my.names[ ,1:t])), ignore.case = T, value = T, invert = F),
               "\n\nYou can also use print.names.under.taxon() to find spellings under a known spelling, for ex., taxon = \"Bacteria\"",
               "will return all the phylum names under domain Bacteria.\nNote taxon names are case sensitive.\n") 
           )
  }
  
  cat(taxon, "is a", colnames(my.names)[taxa.level], "which is taxon level", taxa.level, "\n")
  print(unique(my.names[index, 1:taxa.level]))
  return(taxa.level)
}

#' Print daughter taxon names
#' 
#' For the provided taxon name, print the names of its constituent taxa. This can be used to step down a taxonomy and find the correct spellings of daughter taxa.
#' @param my.list Either a "tax.list" or a "flat.list"
#' @param taxon A correctly spelled and case sensitive taxon name.
#' @param lower.lvl The numeric lowest taxon level names to print. By default, one level below the provided taxon name is printed.
#' @param upper.lvl The numeric upper taxon level names to print. By default no parent taxa are displayed.
#' @param show.lvls FALSE or numeric vector. If a vector taxon levels is provided, it will print the taxon at those levels.
#' @param sort.by How to order the daughter taxa. Either a function such as max or mean, or the character strings "presence" or "alphabet." By default results are shown alphabetically.
#' @export
print.names.under.taxon <- function(my.list, taxon, lower.lvl = "print minimum", upper.lvl = F, show.lvls = F, sort.by = "alphabet"){
  # my.list can be tax.list or flat.list
  
  # if lower.lvl = 8 it prints the full names of all OTUs belonging to that taxon
  # if lower.lvl not specified, it prints 1 level below the input taxon name
  # if upper.lvl = 1 it prints the full taxonomy, including domain. (upper.lvl = 2, starts w phylum)
  # if upper.lvl not specified, it doesn't print any upper level names
  # if specify show.lvls as a vector, it shows those specific levels and overrides lower/upper stuff

  if (!is.matrix(my.list$names)){
    my.names <- my.list$names$seqID
  }else{
    my.names <- my.list$names
  }
  
  t.lvl <- get.taxon.level(my.list = my.list, taxon = taxon)
  
  if (lower.lvl == "print minimum"){
    lower.lvl <- t.lvl + 1
  }
  if (upper.lvl == FALSE){
    upper.lvl = t.lvl + 1
  }
  if (show.lvls[1] == FALSE){
    lvls <- upper.lvl:lower.lvl
  }else{
    lvls <- show.lvls
  }
  
  index <- which(my.names[ ,t.lvl] == taxon)
  
  if (is.function(sort.by)){
    my.list <- flatten.to.single.level(my.list = my.list, finest.level = max(lvls))
    index <- get.taxon.indexes(my.list = my.list, taxa = taxon)
    my.list <- subset.by.taxa(my.list = my.list, keep.index = index, verbose = F, renormalize = F)
    my.list <- sort.by.abundance(my.list = my.list, sort.by = sort.by)
    my.names <- my.list$names[ ,lvls]
  }else if(sort.by == "presence"){
    my.list <- flatten.to.single.level(my.list = my.list, finest.level = max(lvls))
    index <- get.taxon.indexes(my.list = my.list, taxa = taxon)
    my.list <- subset.by.taxa(my.list = my.list, keep.index = index, verbose = F, renormalize = F)
    my.list <- sort.by.abundance(my.list = my.list, sort.by = sort.by)
    my.names <- my.list$names[ ,lvls]
  }else{ # alphabetical
    my.names <- unique(my.names[index, lvls])
  }
  
  return(my.names)
}

#' Get taxon indexes
#' 
#' Return a list of row indexes that correspond to a given taxon or vector of taxa.
#' @param my.list A "flat.list" that includes these names: "names", "av", "sd", "bq", "pd", "br", each of which is an individual matrix.
#' @param taxa A character vector of taxa names. Spelling is case sensitive and must be exact.
#' @return A numeric vector of row indexes corresponding to the taxa provided. This can be input into the subset.by.taxa function.
#' @export
get.taxon.indexes <- function(my.list, taxa){
  # taxon is a character vector of taxa names 
  
  my.names <- my.list$names

  otu.index <- NULL
  for (n in 1:length(taxa)){
    t.lvl <- get.taxon.level(taxon = taxa[n], my.list = my.list)
    index <- which(my.names[ ,t.lvl] == taxa[n])
    if (length(index) < 1){
      t.lvl <- 8
      index <- which(my.names[ ,t.lvl] == taxa[n]) # check if it's an OTU name if it's not a taxa name
    }
    otu.index <- c(otu.index, index)
    cat("keeping", length(index), "OTUs belonging to", colnames(my.names)[t.lvl], taxa[n], "\n\n")
  }
  
  
  keep.index <- unique(otu.index)
  cat("Returning the row indexes for", length(keep.index), "total OTUs that belong to taxa you want to keep.\n")
  if(length(otu.index) > length(keep.index)){
    cat("Note there was overlap in the names you chose to keep.\n")
  }
  
  return(keep.index)
}

#' Normalize abundances by day total or by taxon max values.
#' 
#' @param my.list Either a "tax.list" or a "flat.list"
#' @param by.day Logical. Normalize abundances by sample date- each column divided by its sum. Columns will sum to 100%.
#' @param by.taxon.max Logical. Normalize abundances by taxon max value- each row divided by its max. Rows will have a max value of 100%.
#' @param verbose Logical. Print messages about what's being done.
#' @return Either a "tax.list" or a "flat.list" that has been normalized.
#' After subsetting taxa
#' @export
renormalize <- function(my.list, by.day = TRUE, by.taxon.max = FALSE, verbose = TRUE){
  
  renorm.by.day <- function(my.mat, col.tots){
    my.mat <- t(my.mat)
    my.mat <- my.mat / col.tots * 100
    my.mat <- t(my.mat)
    return(my.mat)
  }
  
  renorm.by.taxon <- function(my.mat, row.tots){
    my.mat <- my.mat / row.tots * 100
    return(my.mat)
  }
  
  # my.list can be tax.list or flat.list or y.vals (flat.list + extra y.max element)
  if ( !is.matrix(my.list$bq) ){ # tax.list
    names(my.list)
    for (t in 1:length(my.list$bq)){
      if (by.day){
        day.tots <- colSums(my.list$av[[t]])
        index.none <- which(day.tots == 0)
        my.list$av[[t]] <- renorm.by.day(my.mat = my.list$av[[t]], col.tots = day.tots)
        my.list$sd[[t]] <- renorm.by.day(my.mat = my.list$sd[[t]], col.tots = day.tots)
        my.list$av[[t]][ ,index.none] <- 0 # if no observations, normalized abundance = 0 not NaN
      }
      if (by.taxon.max){
        taxon.max <- rowSums(my.list$av[[t]])
        index.none <- which(taxon.max == 0)
        my.list$av[[t]] <- renorm.by.taxon(my.mat = my.list$av[[t]], row.tots = taxon.max)
        my.list$sd[[t]] <- renorm.by.taxon(my.mat = my.list$sd[[t]], row.tots = taxon.max)
        my.list$av[[t]][index.none, ] <- 0 # if no observations, normalized abundance = 0 not NaN
        if (verbose){
          cat("\nRenormalizing the average and standard deviation abundance values so max observed abundance is 100% for each taxon.\n")
          cat("The old (not re-normalized) taxon maximums at level ",t," were:\n")
          print(summary(taxon.max))
        }
      }
    }
  }else{ # flat.list or y.vals
    if (by.day){
      day.tots <- colSums(my.list$av)
      index.none <- which(day.tots == 0)
      my.list$av <- renorm.by.day(my.mat = my.list$av, col.tots = day.tots)
      my.list$sd <- renorm.by.day(my.mat = my.list$sd, col.tots = day.tots)
      my.list$av[ ,index.none] <- 0 # if no observations, normalized abundance = 0 not NaN
      if (any(duplicated(c(names(my.list), "filled.sd")))){ # this element exists in y.vals data structure
        my.list$filled.sd <- renorm.by.day(my.mat = my.list$filled.sd, col.tots = day.tots) # note any that are abund zero are no longer filled
      }
    }
    if (by.taxon.max){
      taxon.max <- rowSums(my.list$av)
      index.none <- which(taxon.max == 0)
      my.list$av <- renorm.by.taxon(my.mat = my.list$av, row.tots = taxon.max)
      my.list$sd <- renorm.by.taxon(my.mat = my.list$sd, row.tots = taxon.max)
      my.list$av[index.none, ] <- 0 # if no observations, normalized abundance = 0 not NaN
      if (any(duplicated(c(names(my.list), "filled.sd")))){ # this element exists in y.vals data structure
        my.list$filled.sd <- renorm.by.taxon(my.mat = my.list$filled.sd, row.tots = taxon.max) # note any that are abund zero are no longer filled
      }
      if (verbose){
        cat("Renormalizing the average and standard deviation abundance values so max observed abundance is 100% for each taxon.\n")
        cat("The old (not re-normalized) taxon maximums were:\n")
        print(summary(taxon.max))
      }
    }
  }
  
  if (verbose){
    if (by.day){
      cat("Renormalizing the average and standard deviation abundance values so total OTU abundance is 100% on each day.\n")
      cat("The old (not re-normalized) day totals were:\n")
      print(summary(day.tots))
    }
  }
  
  return(my.list)
}

#' Subset by taxa
#' 
#' Subset a flat.list to include only the provided row indexes. Indexes can be generated by the get.taxon.indexes function. After subsetting, a tax.list can be regenerated with the group.by.tax function.
#' @param my.list A "flat.list" that includes these names: "names", "av", "sd", "bq", "pd", "br", each of which is an individual matrix.
#' @param keep.index A numeric vector of row indexes to keep.
#' @param renormalize Logical. Should samples be re-normalized to that abundances on a given day sum to 100%.
#' @param verbose Logical. 
#' @export
subset.by.taxa <- function(my.list, keep.index, renormalize = F, verbose = T){
  
  # functions ----
  
  subset.to.row.indexes <- function(index.keep, my.list){
    # my.list can be tax.list or flat.list or y.vals (flat.list + extra y.max element)
    if ( !is.matrix(my.list$bq) ){
      names(my.list)
      for (t in 1:length(my.list$bq)){
        my.list$av[[t]] <- my.list$av[[t]][keep.index, ,drop = F]
        my.list$sd[[t]] <- my.list$sd[[t]][keep.index, ,drop = F]
        my.list$bq[[t]] <- my.list$bq[[t]][keep.index, ,drop = F]
        my.list$pd[[t]] <- my.list$pd[[t]][keep.index, ,drop = F]
        my.list$br[[t]] <- my.list$br[[t]][keep.index, ,drop = F]
        my.list$names[[t]] <- my.list$names[[t]][keep.index, ,drop = F]
      }
    }else{
      my.list$av <- my.list$av[keep.index, ,drop = F]
      my.list$sd <- my.list$sd[keep.index, ,drop = F]
      my.list$bq <- my.list$bq[keep.index, ,drop = F]
      my.list$pd <- my.list$pd[keep.index, ,drop = F]
      my.list$br <- my.list$br[keep.index, ,drop = F]
      my.list$names <- my.list$names[keep.index, ,drop = F]
      if (any(duplicated(c(names(my.list), "filled.sd")))){ # this element exists in y.vals data structure
        my.list$filled.sd <- my.list$filled.sd[keep.index, ,drop = F]
        my.list$col.key <- my.list$col.key[keep.index, ,drop = F]
      }
    }
    return(my.list)
  }
  
  renorm <- function(my.mat, col.tots){
    my.mat <- t(my.mat)
    my.mat <- my.mat / col.tots * 100
    my.mat <- t(my.mat)
    return(my.mat)
  }
  
  renormalize.abundances <- function(my.list, renormalize, message){
    
    if (message){
      new.tots <- colSums(my.list$av)
      cat("You are subsetting your samples to this percent of the orig. reads:\n")
      print(summary(new.tots))
      cat("\n")
    }else if(renormalize){
      new.tots <- colSums(my.list$av)
    }
    
    if (renormalize){
      cat("Renormalizing the average and standard deviation abundance values so total OTU abundance is 100% on each day.\n")
      my.list$av <- renorm(my.mat = my.list$av, col.tots = new.tots)
      my.list$sd <- renorm(my.mat = my.list$sd, col.tots = new.tots)
    }
    
    return(my.list)
  }
  
  # actions ----
  
  my.list <- subset.to.row.indexes(index.keep = keep.index, my.list = my.list)
  
  my.list <- renormalize.abundances(my.list = my.list, renormalize = renormalize, message = verbose)
  
  return(my.list)
}

#' Make abundances below the limit of quantification equal zero
#'
#' If a limit of quantification definition is chosen, this will make all values below it equal to zero. The LOQ can be defined as below resolution (br): below to abundance of a single read in the shallowest sample, partially detected (pd): not detected in all replicates, or below quantification (bq): either br or pd.
#' @param my.list Either a "tax.list" or a "flat.list".
#' @param LOQ.def Either "bq", "br", or "pd"
#' @param renorm Logical. Should samples dates be renomralized to 100% total abundance.
#' @param verbose Logical.
#' @export
make.zero.below.LOQ <- function(my.list, LOQ.def = "bq", renorm = TRUE, verbose = TRUE){
  # my.list can be tax.list or flat.list or y.vals (flat.list + extra y.max element)
  if ( !is.matrix(my.list$bq) ){
    names(my.list)
    for (t in 1:length(my.list$bq)){
      my.list$av[[t]] <- my.list$av[[t]] * !my.list[[LOQ.def]][[t]]
      my.list$sd[[t]] <- my.list$sd[[t]] * !my.list[[LOQ.def]][[t]]
    }
  }else{
    my.list$av <- my.list$av * !my.list[[LOQ.def]]
    my.list$sd <- my.list$sd * !my.list[[LOQ.def]]
    
    if (any(duplicated(c(names(my.list), "filled.sd")))){ # this element exists in y.vals data structure
      my.list$filled.sd <- my.list$filled.sd * !my.list[[LOQ.def]]
    }
  }
  
  if (renorm){
    my.list <- renormalize(my.list = my.list, by.day = TRUE, by.taxon.max = FALSE, verbose = verbose)
  }
  
  return(my.list)
}


