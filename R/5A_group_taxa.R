# RRR 10-14-20

# ---- sourceable grouping function! ----

#' Make an empty list structure.
#' 
#' Make a where each element is NULL but has the provided names.
#' @param ListNames A vector of names to use
#' @return A list where each element is NULL but is named with the ListNames
#' @export
make.empty.list.structure <- function(ListNames){
  # the ListNames can be something like c("OTU", "kingdom","phylum","class","order","family/lineage","genus/clade","species/tribe")
  empty.list <- list(NULL)
  for (e in 1:length(ListNames)){
    empty.list[[e]] <- 0
    names(empty.list)[e] <- ListNames[e]
  }
  return(empty.list)
}

#' Group by taxonomy
#' 
#' Group a "flat list" structure by taxonomy, so that it becomes a nested list with each taxon level included separately. In addition to summing daughter taxa abundances, this function propagates error for the sd and calculates quantification limits.
#' @param my.list A "flat.list" that includes these names: "names", "av", "sd", "bq", "pd", "br", each of which is an individual matrix.
#' @return A "tax.list" that has these names: "names", "av", "sd", "bq", "pd", "br", each of which is a nested list that has a matrix for each taxonomy level.
#' @export
group.by.tax <- function(my.list){
  # my.list is either the combo.list or the subset.list
  # should include these names, which contain individual matrices
  #   "names"    "av"    "sd"    "bq"    "pd"    "br"
  
  # ---- functions ----
  
  propagate.error.sums <- function(error.vector){
    # note: there's a package called propagate that can figure out more complex formulas. seems easier to do by hand since formulas simple.
    sx <- sqrt( sum( error.vector ^ 2 ) )
    return(sx)
  }
  
  group <- function(my.names, my.vals, my.fun, my.levels){
    # the aggregation list contains all of the taxa names (not unique) for each level 
    # aggregate will combine rows of abundance that have the same taxa name
    
    aggregation.list <- make.empty.list.structure(ListNames = my.levels)
    for(t in 1:length(aggregation.list)){
      aggregation.list[[t]] <- my.names[ ,t]
    }
    
    grouped.list <- make.empty.list.structure(ListNames = my.levels)
    for (t in 1:length(grouped.list)){
      grouped.list[[t]] <- aggregate(x = my.vals, by = aggregation.list[1:t], FUN = my.fun)
    }
    
    return(grouped.list)
  }
  
  clean.up.list <- function(grouped.list, lowest.level.vals, lowest.level.names, my.levels, pull.out.names = FALSE){
    
    # put grouped taxa rows in alphabetical order
    for(t in 1:length(grouped.list)){
      index <- do.call(what = order, args = grouped.list[[t]][ ,1:t, drop = F])
      grouped.list[[t]] <- grouped.list[[t]][index, ]
    }
    
    # put the OTU-level that wasn't aggregated in alphabetical order too
    index <- do.call(what = order, args = as.data.frame(lowest.level.names))
    lowest.level.names <- lowest.level.names[index, ,drop = FALSE]
    lowest.level.vals <- lowest.level.vals[index, ,drop = FALSE]
    # cat(all.equal(row.names(lowest.level.vals), lowest.level.names[,ncol(lowest.level.names)]))
    
    # keep only names or only values, convert to matrix
    if (pull.out.names == FALSE){
      for(t in 1:length(grouped.list)){
        save.row.names <- grouped.list[[t]][ ,t]
        grouped.list[[t]] <- grouped.list[[t]][ ,-(1:t), drop = F]
        grouped.list[[t]] <- as.matrix(grouped.list[[t]])
        row.names(grouped.list[[t]]) <- save.row.names
      }
      lowest.level <- lowest.level.vals
      
    }else if(pull.out.names == TRUE){
      for(t in 1:length(grouped.list)){
        grouped.list[[t]] <- grouped.list[[t]][ ,(1:t), drop = F]
        grouped.list[[t]] <- as.matrix(grouped.list[[t]])
        row.names(grouped.list[[t]]) <- NULL
      }
      lowest.level <- lowest.level.names
    }
    
    # add in the un-aggregated OTU level (already matrix)
    next.level <- length(grouped.list) + 1
    grouped.list[[next.level]] <- lowest.level
    names(grouped.list)[next.level] <- my.levels[next.level]
    
    return(grouped.list)
  }
  
  check.name.order <- function(){ # lazy calls to global env
    for (t in 1:length(names(names.list))){
      cat(t,"-", names(names.list)[t], ": ", 
          all.equal(as.vector(names.list[[t]][ ,t]), row.names(av.list[[t]])),
          all.equal(as.vector(names.list[[t]][ ,t]), row.names(sd.list[[t]])),
          all.equal(as.vector(names.list[[t]][ ,t]), row.names(bq.list[[t]])),
          all.equal(as.vector(names.list[[t]][ ,t]), row.names(pd.list[[t]])),
          all.equal(as.vector(names.list[[t]][ ,t]), row.names(br.list[[t]])),
          "\n")
    }
  }
  
  
  # ---- group by taxa ----
  
  taxa.levels <- c("Kingdom","Phylum","Class","Order","Family/Lineage","Genus/Clade","Species/Tribe","seqID")
  
  if (ncol(my.list$names) < 8){
    agg.levels <- taxa.levels[1:ncol(my.list$names)]
    included.levels <- agg.levels
    agg.levels <- agg.levels[-length(agg.levels)]
  }else{
    included.levels <- taxa.levels
    agg.levels <- taxa.levels[-8] # it's WAY faster to not run aggregate the OTU level, which is already unique anyway.  
  }
 
  
  cat("\nCombining mean abundances by summing taxa.\n")
  av.list <- group(my.names = my.list$names, my.vals = my.list$av, my.fun = sum, my.levels = agg.levels)
  
  cat("Calculating standard deviation by propagating error of sums.\n")
  sd.list <- group(my.names = my.list$names, my.vals = my.list$sd, my.fun = propagate.error.sums, my.levels = agg.levels)
  
  cat("Classifying as below quantification if not all constituent taxa were above quantification.\n")
  bq.list <- group(my.names = my.list$names, my.vals = my.list$bq, my.fun = all, my.levels = agg.levels)
  pd.list <- group(my.names = my.list$names, my.vals = my.list$pd, my.fun = all, my.levels = agg.levels)
  br.list <- group(my.names = my.list$names, my.vals = my.list$br, my.fun = all, my.levels = agg.levels)
  
  cat("Putting all the rows in alphabetical order by taxonomy.\n")
  names.list <- clean.up.list(grouped.list = bq.list, lowest.level.vals = my.list$bq, lowest.level.names = my.list$names, my.levels = included.levels, 
                              pull.out.names = TRUE)
  av.list <- clean.up.list(grouped.list = av.list, lowest.level.vals = my.list$av, lowest.level.names = my.list$names, my.levels = included.levels)
  sd.list <- clean.up.list(grouped.list = sd.list, lowest.level.vals = my.list$sd, lowest.level.names = my.list$names, my.levels = included.levels)
  bq.list <- clean.up.list(grouped.list = bq.list, lowest.level.vals = my.list$bq, lowest.level.names = my.list$names, my.levels = included.levels)
  pd.list <- clean.up.list(grouped.list = pd.list, lowest.level.vals = my.list$pd, lowest.level.names = my.list$names, my.levels = included.levels)
  br.list <- clean.up.list(grouped.list = br.list, lowest.level.vals = my.list$br, lowest.level.names = my.list$names, my.levels = included.levels)
  
  cat("\ncheck row orders match:\n")
  cat(check.name.order())
  cat("check if abundances normalized by sample:\n")
  for (t in 1:length(included.levels)){
    cat(summary(colSums(av.list[[t]][ ,-(1:t), drop = F])),"\n") # this should be true from previous scripts, shouldn't need to re-nomalize here!
  }
  
  # the pvals info only applies to the OTU-level, once combine OTUs you're combining names with diff pvalues
  # no need to hold onto this anymore, just go back to combo.list for it
  
  # the resolution only applies to the individual replicates. Once replicates are combined there is no longer any way to back-calculate to absolute reads.
  # no need to hold onto this anymore, just go back to combo.list for a summary or taxass.samples for direct comparison
  # the br matrix and the pd matrix carry this information through
  
  tax.list <- list("names" = names.list,
                   "av" = av.list,
                   "sd" = sd.list,
                   "bq" = bq.list,
                   "pd" = pd.list,
                   "br" = br.list)
  
  return(tax.list)
}


# ---- end ----

# # ---- Checked this example's math by hand ----
# 
# test.matrix <- matrix(nrow = 10, ncol = 4, sample(x = 0:5, size = 40, replace = T))
# 
# test.names <- matrix(nrow = 10, ncol = 2, data = c("A","A","A","A","B","B","B","B","B","B", letters[1:10]))
# 
# test.list <- group.sequences(Taxonomy = test.names, Abunds = test.matrix, ListNames = c("tribe","otu"))
# 
# test.sds <- prop.error.group.sequences(Taxonomy = test.names, Abunds = test.matrix, ListNames = c("tribe","otu"))
# 
# test.list
# test.sds

# # ---- used to test finding zeros and NAs ----
# 
# # test matrix creation:
# a.matrix <- matrix(data = 5, nrow = 4, ncol = 6)
# a.matrix[c(2,4,6,10,15,16,20,24)] <- 0
# b.matrix <- matrix(data = 6, nrow = 4, ncol = 6)
# b.matrix[c(1:3,5,7,8,14,16,18,21)] <- 0
# b.matrix[ ,2] <- NA
