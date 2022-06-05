# RRR 
# Combine bological replicates into mean, standard deviation, and below quantification tables

# ---- Set-Up ----

current.data.folder <- "data/limony-IRD/2021-08-25_processing"

input.samples.list <- "3A_taxass_samples_6000_pn0.05_t1.rds"

input.samples.key <- "2C_key_samples_6000.rds"

output.combo.list <- sub(pattern = "taxass", replacement = "combo", x = input.samples.list)
output.combo.list <- sub(pattern = "2C", replacement = "4A", x = output.combo.list)
output.combo.list <- sub(pattern = "3A", replacement = "4A", x = output.combo.list)
output.combo.list

library("lubridate")


# ---- Functions ----

top <- function(X, r = 5, c = 10, right = F, bottom = F){
  # not table
  if(!is.data.frame(X) & !is.matrix(X)){
    
    # vector
    if (!is.list(X)){ 
      if (r >= length(X)){
        r <- length(X)
      }else{
        cat("length = ", length(X), "\n")
      }
      if(bottom | right){
        rows <- (length(X) - r + 1):length(X)
      }else{
        rows <- 1:r
      }
      return(X[rows])
    }
    
    # list
    cat("last element of list with length = ", length(X), "\n")
    X <- X[[length(X)]]
  }
  
  # matrix or data frame
  if(r > nrow(X)){
    r <- nrow(X)
  }
  if(c > ncol(X)){
    c <- ncol(X)
  }
  
  if(right){
    cols <- (ncol(X) - r + 1):ncol(X)
  }else{
    cols <- 1:c
  }
  
  if(bottom){
    rows <- (nrow(X) - r + 1):nrow(X)
  }else{
    rows <- 1:r
  }
  if(nrow(X) > r){
    cat("nrow = ", nrow(X), "\n")
  }
  if(ncol(X) > c){
    cat("ncol = ", ncol(X), "\n")
  }
  return(X[rows,cols])
}

make.empty.list.structure <- function(ListNames){
  # the ListNames can be something like c("OTU", "kingdom","phylum","class","order","family/lineage","genus/clade","species/tribe")
  empty.list <- list(NULL)
  for (e in 1:length(ListNames)){
    empty.list[[e]] <- 0
    names(empty.list)[e] <- ListNames[e]
  }
  return(empty.list)
}

get.all.duplicated.indexes <- function(my.vector){
  index1 <- which(duplicated(x = my.vector, incomparables = NA))
  index2 <- which(duplicated( x = c(my.vector[index1], my.vector), incomparables = NA ))
  index2 <- index2 - length(index1)
  return(index2)
}

convert.sample.names.to.dates <- function(sample.names){
  sample.names <- sub(pattern = "\\..*$", replacement = "", x = sample.names)
  sample.names <- parse_date_time(x = sample.names, orders = "dmy", tz = "Etc/GMT-5")
  return(sample.names)
}

get.rep.index.table <- function(abunds, key){
  new.names <- data.frame("limony.names" = colnames(abunds), "order" = 1:ncol(abunds))
  key <- data.frame("limony.names" = key$limony.names, "in.R.colnames" = key$in.R.colnames, "Lim.Rep" = key$Lim.Rep)
  new.names <- merge(x = new.names, y = key, by = "limony.names", all = TRUE)
  index <- order(new.names$order)
  new.names <- new.names[index, ]
  return(new.names)
}

separate.by.replicate <- function(taxass, key){
  my.reps <- unique(key$Lim.Rep)
  my.reps <- sort(my.reps)
  abunds.split <- make.empty.list.structure(ListNames = my.reps)
  res.split <- make.empty.list.structure(ListNames = my.reps)
  num.samples <- length(unique(key$in.R.colnames))
  for (r in 1:length(abunds.split)){
    index <- which(key$Lim.Rep == names(abunds.split)[r])
    abunds.split[[r]] <- taxass$abunds[ ,index] 
    res.split[[r]] <- taxass$res[index]
  }
  
  return(list("abund" = abunds.split, "res" = res.split))
}

rename.sample.columns <- function(reps, key){
  index <- which(colnames(key) == "order")
  key <- key[ ,-index]
  
  for (r in 1:length(reps$abund)){
    new.names <- data.frame("limony.names" = colnames(reps$abund[[r]]), "order" = 1:ncol(reps$abund[[r]]))
    new.names <- merge(x = new.names, y = key, by = "limony.names", all.x = TRUE)
    index <- order(new.names$order)
    new.names <- new.names[index, ]
    
    cat(all.equal(colnames(reps$abund[[r]]), new.names$limony.names),
        all.equal(names(reps$res[[r]]), new.names$limony.names))
    colnames(reps$abund[[r]]) <- new.names$in.R.colnames
    names(reps$res[[r]]) <- new.names$in.R.colnames
  }
  return(reps)
}

pool.sequencing.replicates <- function(reps){
  for (r in 1:length(reps$abund)){
    index <- duplicated(colnames(reps$abund[[r]]))
    seq.reps <- unique(colnames(reps$abund[[r]])[index])
    
    if (length(seq.reps) < 1){
      cat("(There are no sequencing reps in Bio Rep", names(reps$abund)[r], ")\n")
      next
    }
    
    for (s in 1:length(seq.reps)){
      index <- which(colnames(reps$abund[[r]]) == seq.reps[s])
      cat("Combining (pooling)", length(index), "sequencing replicates for", seq.reps[s], "(Bio Rep", names(reps$abund)[r], ")\n")
      
      pooled.abund <- reps$abund[[r]][ ,index]
      pooled.abund <- rowSums(pooled.abund)
      pooled.abund <- pooled.abund / sum(pooled.abund) * 100 # re-close to 100%
      
      pooled.res <- reps$res[[r]][index]
      pooled.res <- sum(pooled.res)
      
      first.rep <- index[1]
      other.reps <- index[-1]
      
      reps$abund[[r]][ ,first.rep] <- pooled.abund
      reps$abund[[r]] <- reps$abund[[r]][ ,-other.reps]
      
      reps$res[[r]][first.rep] <- pooled.res
      reps$res[[r]] <- reps$res[[r]][-other.reps]
    }
  }
  return(reps)
}

add.missing.sample.columns <- function(reps, key){
  
  # first add columns
  for (r in 1:length(reps$abund)){
    missing.names <- setdiff(x = key$in.R.colnames, y = colnames(reps$abund[[r]]))
    if(length(missing.names) < 1){next}
    
    temp <- matrix(data = NA, nrow = nrow(reps$abund[[r]]), ncol = length(missing.names))
    colnames(temp) <- missing.names
    temp <- as.data.frame(temp)
    reps$abund[[r]] <- cbind(reps$abund[[r]], temp)
    reps$abund[[r]] <- as.matrix(reps$abund[[r]])
    
    temp <- rep.int(x = NA, times = length(missing.names))
    names(temp) <- missing.names
    reps$res[[r]] <- c(reps$res[[r]], temp)
  }
  
  # then order columns
  for (r in 1:length(reps$abund)){
    sample.dates <- convert.sample.names.to.dates(sample.names = colnames(reps$abund[[r]]))
    index <- order(sample.dates)
    reps$abund[[r]] <- reps$abund[[r]][ ,index]
    reps$res[[r]] <- reps$res[[r]][index]
    cat(all.equal(colnames(reps$abund[[r]]), names(reps$res[[r]])))
  }
  cat("check")
  cat(all.equal(colnames(reps$abund$A), colnames(reps$abund$B)))
  cat(all.equal(colnames(reps$abund$A), colnames(reps$abund$C)))
  return(reps)
}

combine.rep.resolutions <- function(res){
  res.tab <- matrix(data = NA, ncol = length(res[[1]]), nrow = length(res))
  colnames(res.tab) <- names(res[[1]])
  row.names(res.tab) <- names(res)
  for (r in 1:length(res)){
    res.tab[r, ] <- res[[r]]
  }
  return(res.tab)
}

calc.res.as.perc.abund <- function(res){
  # the percent abundance version of res will be useful later, because that can be carried through when subsetting
  # while reads/run cannot be recalculated after re-normalizing b/c don't know respective contribution to average of each replicate
  # the abundance of 1 read after subsetting can be approximately tracked, by normalizing to the same new total abundance
  # could similarly take a percent of read counts to adjust res, but that gives a false certainty of reads/run that are still included
  # note: this also has an extra row- the daily worst
  
  res <- 1 / res * 100
  worst <- apply(X = res, MARGIN = 2, FUN = min, na.rm = T)
  res <- rbind(res, worst)
  
  return(res)
}

calc.mean.of.reps <- function(reps, use.N.minus.1 = FALSE){
  # reps is reps$abund not rep.list
  
  for (r in 1:length(reps)){
    x <- reps[[r]]
    na.index <- is.na(x)
    x.s <- x
    x.s[na.index] <- 0
    x.d <- x
    x.d[na.index] <- 0
    x.d[!na.index] <- 1
    reps[[r]] <- list("sum" = x.s, "div" = x.d)
  }
  
  cum.sum <- matrix(data = 0, nrow = nrow(reps[[1]][[1]]), ncol = ncol(reps[[1]][[1]]))
  for (r in 1:length(reps)){
    cum.sum <- cum.sum + reps[[r]]$sum
  }
  
  cum.denom <- matrix(data = 0, nrow = nrow(reps[[1]][[1]]), ncol = ncol(reps[[1]][[1]]))
  for (r in 1:length(reps)){
    cum.denom <- cum.denom + reps[[r]]$div
  }
  
  if (use.N.minus.1){
    cum.denom <- cum.denom - 1
  }
  
  av <- cum.sum / cum.denom
  return(av)
}

calc.sd.of.reps <- function(reps, av){
  # sd is the MEAN of the squared differenced. 
  # so just change rep values to be (x-m)^2 and then reuse my mean function that ignores NAs
  # BUT for sd you use N-1 as the denominator when taking the mean, instead of N
  
  for (r in 1:length(reps)){
    x <- reps[[r]]
    x <- (x - av)^2
    reps[[r]] <- x
  }
  
  my.sd <- calc.mean.of.reps(reps = reps, use.N.minus.1 = TRUE)
  my.sd <- sqrt(my.sd)
  
  # when only have 1 rep, SD is NaN. change it to NA
  index <- is.nan(my.sd)
  my.sd[index] <- NA
  
  return(my.sd)
}

calc.if.partial.detection <- function(reps){
  # If present in some but not all reps --> TRUE
  # If present in all reps --> FALSE
  # If zero in all reps --> NA
  # If only 1 rep --> NA
  
  for (r in 1:length(reps)){
    x <- reps[[r]]
    na.index <- is.na(x)
    abs.index <- x == 0
    
    x.r <- x  # does replicate exist?
    x.r[na.index] <- 0
    x.r[!na.index] <- 1
    
    x.a <- x # was OTU absent
    x.a[abs.index] <- TRUE
    x.a[!abs.index] <- FALSE
    x.a[na.index] <- FALSE 
    
    reps[[r]] <- list("is.absent" = x.a, "has.replicate" = x.r)
  }
  
  num.reps.exist <- matrix(data = 0, nrow = nrow(reps[[1]][[1]]), ncol = ncol(reps[[1]][[1]]))
  for (r in 1:length(reps)){
    num.reps.exist <- num.reps.exist + reps[[r]]$has.replicate
  }
  num.reps.exist
  
  is.it.ever.absent <- matrix(data = FALSE, nrow = nrow(reps[[1]][[1]]), ncol = ncol(reps[[1]][[1]]))
  for (r in 1:length(reps)){
    is.it.ever.absent <- is.it.ever.absent + reps[[r]]$is.absent
  }
  is.it.ever.absent
  
  is.it.always.absent <- num.reps.exist == is.it.ever.absent
  
  do.reps.exist <- num.reps.exist > 1
  
  is.ud <- is.it.ever.absent > 0
  is.ud[!do.reps.exist] <- NA
  is.ud[is.it.always.absent] <- NA
  return(is.ud)
}

calc.if.below.resolution <- function(reps, res){
  # below resolution if ANY value is below resolution
  # do by each replicate & then combine. if only 1 rep, then br based on that one only
  
  res <- min(res, na.rm = T)
  res.perc <- 1 / res *100
  cat("The worst resolution sample has", res, "reads, which means the resolution limit is 1 /", res, "(*100%) =", res.perc, "% abundance.\n")
  
  for (r in 1:length(reps)){
    x <- reps[[r]]
    
    x.br <- x < res.perc
    index.na <- is.na(x.br)
    x.br[index.na] <- FALSE # when the rep doesn't exist, it's not below resolution
    
    reps[[r]] <- x.br
  }
  
  cum.br <- matrix(data = FALSE, nrow = nrow(reps[[1]]), ncol = ncol(reps[[1]]))
  for (r in 1:length(reps)){
    x <- reps[[r]]
    
    cum.br <- x | cum.br
  }
  
  return(cum.br)
}

calc.below.quantification <- function(pd, br){
  # bq if average value is less that the worst resolution of a single read
  # bq if there were replicates and it was not present in all of them
  # if there were no replicates, it's based on abundance only. there are no NA's in bq
  
  na.index <- is.na(pd)
  pd[na.index] <- FALSE # being NA b/c all 0 is TRUE in br anyway, not bq just because unduplicated
  
  bq <- pd | br
  return(bq)
}


# ---- Go ----

taxass <- readRDS(file = file.path(current.data.folder, input.samples.list))

key <- readRDS(file = file.path(current.data.folder, input.samples.key))

key <- get.rep.index.table(abunds = taxass$abunds, key = key)

rep.list <- separate.by.replicate(taxass = taxass, key = key)

rep.list <- rename.sample.columns(reps = rep.list, key = key)

rep.list <- pool.sequencing.replicates(reps = rep.list)

rep.list <- add.missing.sample.columns(reps = rep.list, key = key)

combo.list <- make.empty.list.structure(ListNames = c("names", "av","sd","bq","pd","br", "res.perc", "res.reads", "pvals")) # "average" (mean prercent), "standard deviation", "below quantification" means either "partial detection" in some but not all of the replicates or "below resolution" of the worst-resolution sample

combo.list$res.reads <- combine.rep.resolutions(res = rep.list$res)

combo.list$res.perc <- calc.res.as.perc.abund(res = combo.list$res.reads)

combo.list$av <- calc.mean.of.reps(reps = rep.list$abund)

combo.list$sd <- calc.sd.of.reps(reps = rep.list$abund, av = combo.list$av)

combo.list$pd <- calc.if.partial.detection(reps = rep.list$abund)

combo.list$br <- calc.if.below.resolution(reps = rep.list$abund, res = combo.list$res.reads)

combo.list$bq <- calc.below.quantification(pd = combo.list$pd, br = combo.list$br)

combo.list$names <- taxass$names

combo.list$pvals <- taxass$pvalues

# ---- Export combined replicates list ----

output.file <- file.path(current.data.folder,output.combo.list)
cat("Making file: ", output.file, "\n")
# saveRDS(object = combo.list, file = output.file)



# # ---- Test calcs ----
# x = matrix(data = c(3,4,0,0,5,1), nrow = 2, ncol = 3, byrow = T)
# y = matrix(data = c(1,6,NA,0,0,NA), nrow = 2, ncol = 3, byrow = T)
# z = matrix(data = c(4,NA,NA,0,NA,NA), nrow = 2, ncol = 3, byrow = T)
# x
# y
# z
# reps = list(A = x, B = y, C = z)
# for (r in 1:length(reps)){
#   colnames(reps[[r]]) <- paste0("date",1:ncol(x))
#   rownames(reps[[r]]) <- paste0("otu", 1:nrow(x))
# }
# reps
# 
# av <- calc.mean.of.reps(reps = reps)
# av
# std <-calc.sd.of.reps(reps = reps, av = av)
# std
# pd <- calc.if.partial.detection(reps = reps)
# pd
# br <- calc.if.below.resolution(av = av, res = 50)
# br
# bq <- calc.below.quantification(pd = pd, br = br)
# bq
