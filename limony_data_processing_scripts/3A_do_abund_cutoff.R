# RRR
# Remove low-abundance OTUs from the dataset
# Define abundance based on number reads (n) and/or times observed (t)

# ---- set-up ----

current.data.folder <- "data/limony-IRD/2021-08-25_processing"

input.ungrouped.list <- "2C_taxass_samples_6000.rds"

min.spike.filter <- .1
min.times.filter <- 1
time.as.percent <- FALSE
abund.as.percent <- TRUE
use.sum.abundance <- FALSE                                

library("magrittr")

# created file has format inputname_totpn#_pt#.rds where tot indicated used sum abundance and p indicates used percent

# ---- define functions ----

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

generate.output.filename <- function(cutoffs, my.input, my.folder){
  if (cutoffs$n.is.perc){pn = "p"}else{pn = ""}
  if (cutoffs$t.is.perc){pt = "p"}else{pt = ""}
  if (cutoffs$n.is.tot){sa = "tot"}else{sa = ""}
  no.rds <- sub(pattern = "\\.rds", replacement = "", x = my.input)
  no.rds <- sub(pattern = "^2C_", replacement = "3A_", x = no.rds)
  
  
  created.file <- paste0(no.rds, "_", sa, pn, "n", cutoffs$min.abund, "_", pt, "t", cutoffs$min.times, ".rds")
  created.file <- file.path(my.folder, created.file)
  
  return(created.file)
}

unnormalize.matrix <- function(abund.matrix, res.vector){
  seqID.names <- row.names(abund.matrix)
  abund.matrix <- t(abund.matrix)
  abund.matrix <- abund.matrix * res.vector / 100
  abund.matrix <- t(abund.matrix)
  abund.matrix <- round(x = abund.matrix, digits = 0)
  abund.matrix <- apply(X = abund.matrix, MARGIN = 2, FUN = as.integer)
  row.names(abund.matrix) <- seqID.names
  # all.equal(colSums(abund.matrix), res.vector) #check
  return(abund.matrix)
}

normalize.matrix <- function(abund.matrix){
  tots <- colSums(abund.matrix)
  abund.matrix <- t(abund.matrix) 
  abund.matrix <- abund.matrix / tots * 100
  abund.matrix <- t(abund.matrix)
  
  check <- unique(colSums(abund.matrix))
  check2 <- NULL
  for (c in check){
    check2 <- c(check2, all.equal(c, 100))
  }
  if (!all(check2)){ 
    return(cat("\n\nPROBLEM WITH NORMALIZING BY SAMPLE READS!\n\n"))
  }
  
  return(abund.matrix)
}

print.directions <- function(){
  # cat("\nabund table has seqID in rows, samples in cols, normalized by sample so each col sums to 100 %.\n\n")
  return(cat("To be kept, a seqID must reach the min.abund in at least min.times of samples.\n",
             "   i.e. seqID >= n for >= t times\n\n",
             
             "If instead you want to filter \"at least min.times\" and \"at least one time at min.abund\" use this function sequentially.\n",
             "   set min spike = 0 to filter just on presence/absence (special case: filters with > 0 instead of >= 0)\n",
             "   set min times = 1 and time.as.perc = FALSE to filter just on max abund\n\n", 
             
             "Examples:\n", 
             "   keep everything: \n      min.spike = 0, min.times = 0\n", 
             "   remove seqs that only occur once ever: \n      min.spike = 2, min.times = 1, time.as.perc = F, abund.as.perc = F, overall.abundance = T\n", 
             "   keep things that occur with at least 0.05 % rel abund: \n      min.spike = .05, min.times = 1, abund.as.perc = T, time.as.percent = F\n", 
             "   keep things that occur with at least 0.05 % rel abund in at least 2 samples: \n      min.spike = .05, min.times = 2, abund.as.perc = T, time.as.percent = F\n", 
             "   keep things that occur with at least 0.05 % rel abund in at least 2% of the samples: \n      min.spike = .05, min.times = 1, abund.as.perc = T, time.as.percent = T\n", 
             "   keep things that contribte to at least 1 % of total abundance: \n      min.spike = 1, abund.as.perc = T, overall.abundance = T\n\n", sep = ""))
}

find.rows.above.abund.cutoff <- function(my.list, cutoffs){
  with(my.list, {
    with(cutoffs, {
      
      if (n.is.perc == FALSE){
        abunds <- unnormalize.matrix(abund.matrix = abunds, res.vector = res)
      }
      
      if (n.is.tot == TRUE){
        abunds <- rowSums(abunds)
        abunds <- as.matrix(x = abunds, ncol = 1)
        if(t.is.perc != FALSE | min.times != 1){
          cat("Can't filter by time if using total abundances. Forcing t.is.perc = F and min.times = 1")
          t.is.perc <- FALSE
          min.times <- 1
        }
      }
      
      if (min.abund == 0){
        spike.filtered <- abunds > min.abund
      }else{
        spike.filtered <- abunds >= min.abund
      }
      
      if (t.is.perc == TRUE){       
        min.times <- min.times * ncol(abunds) / 100 
      }
      
      keep.rows.index <- rowSums(spike.filtered) >= min.times
      return(keep.rows.index)
    })
  })
}

apply.cutoff <- function(my.list, keep.index){
  
  my.list <- within(my.list, {
    abunds <- abunds[keep.index, , drop = F]
    names <- names[keep.index, , drop = F]
    pvalues <- pvalues[keep.index, , drop = F]
  })
  
  my.list$res <- with(my.list, {
    unnormalize.matrix(abund.matrix = abunds, res.vector = res) %>%
      colSums(.)
  })
  
  my.list$abunds <- normalize.matrix(abund.matrix = my.list$abunds)
  
  return(my.list)
}

how.many.removed <- function(full.list, keep.list){ 
  
  old.otus <-  nrow(full.list$abunds)
  new.otus <- nrow(keep.list$abunds)
  total.otus.removed <- old.otus - new.otus
  perc.otus.removed <- total.otus.removed / old.otus * 100
  
  old.reads <- (sum(full.list$res))
  new.reads <- sum(keep.list$res)
  total.reads.removed <- old.reads - new.reads 
  perc.reads.removed <- total.reads.removed / old.reads * 100
  
  sum.table <- data.frame(total.otus.removed, perc.otus.removed, total.reads.removed, perc.reads.removed)
  sum.table <- round(sum.table, digits = 2)
  return(sum.table)
}

check.seqID.orders.match <- function(abund.IDs, tax.IDs){
  if (all.equal(abund.IDs, tax.IDs) == TRUE){
    cat("ok good keep going")
  }else{
    cat("fuck fuck fuck fuck fuck\n")
    cat(all.equal(abund.IDs, tax.IDs))
    cat("\n\n")
    cat("abund.IDs[1]\n", abund.IDs[1])
    cat("\n\ntax.IDs[1]\n", tax.IDs[1])
  }
}


# ---- use functions ---- 

taxass <- readRDS(file = file.path(current.data.folder, input.ungrouped.list))

cutoffs <- data.frame("min.times" = min.times.filter, "t.is.perc" = time.as.percent, "min.abund" = min.spike.filter, "n.is.perc" = abund.as.percent, "n.is.tot" = use.sum.abundance)

created.list.filename <- generate.output.filename(cutoffs = cutoffs, my.input = input.ungrouped.list, my.folder = current.data.folder)
created.stats.filename <- sub(pattern = ".rds", replacement = "_cutoff-summary.csv", x = created.list.filename)

print.directions()

index.above.cutoff <- find.rows.above.abund.cutoff(my.list = taxass, cutoffs = cutoffs)

taxass.keep <- apply.cutoff(my.list = taxass, keep.index = index.above.cutoff)

# taxass.toss <- apply.cutoff(my.list = taxass, keep.index = -index.above.cutoff) # If care to look at what is below the cutoff

cutoff.stats <- how.many.removed(full.list = taxass, keep.list = taxass.keep)

cutoff.stats <- cbind(cutoffs, cutoff.stats)
cutoff.stats

# ---- Save files ----

cat("making file: ", created.list.filename, "\n")
# saveRDS(object = taxass.keep, file = created.list.filename)

cat("made file: ", created.stats.filename, "\n")
# write.csv(x = cutoff.stats, file = created.stats.filename, quote = F, row.names = F)


# end

