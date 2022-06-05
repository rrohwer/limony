# RRR 7/12/20
# 1.
# remove the standards (b/c many blanks "fail" sequencing)
# remove samples below sequencing read cutoff (samples that "failed" sequencing)
# left with only the lake.samples that will use in the dataset.
# 2. 
# edit sample names to add a pf after pre-filtered sampled (to distinguish from same-named ww samples)
# edit biological, filter, and extraction replicates to have consistent nomenclature
# assign "limony" rep for how to handle them in this dataset for comparison purposes
# 3.
# Save a fully expanded abundance table, with rr's as column names but with rest of info added to the key.


# ---- set-up ----

CHOSEN.CUTOFF <- 6000

library(lubridate)

current.r.data.folder <- "data/limony-IRD/2021-08-25_processing"
extractions.data.entry.folder <- "data/metadata"

input.formatted.taxass <- "1B_taxass.rds"
input.data.entry.extractions <- "2020-07-10_data_entry_limony_and_TYMEFLIES.rds"

created.samples.list <- "2C_taxass_samples_6000.rds"
created.samples.key <- "2C_key_samples_6000.rds"


# ---- functions ----

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

get.all.duplicated.indexes <- function(my.vector){
  index1 <- which(duplicated(x = my.vector, incomparables = NA))
  index2 <- which(duplicated( x = c(my.vector[index1], my.vector), incomparables = NA ))
  index2 <- index2 - length(index1)
  return(index2)
}

get.key.samples <- function(key){
  index <- which(key$sample.type == "lake.sample" & !is.na(key$limony.names))
  key <- key[index, ]
  rownames(key) <- 1:nrow(key)
  return(key)
}

get.taxass.samples <- function(taxass, key){
  index <- is.element(el = names(taxass$res), set = key$limony.names)
  taxass$res <- taxass$res[index]
  taxass$abunds <- taxass$abunds[ ,index]
  
  cat(all.equal(names(taxass$res), colnames(taxass$abunds)))
  return(taxass)
}

get.taxass.passing <- function(taxass, cutoff){
  index <- taxass$res >= cutoff
  taxass$res <- taxass$res[index]
  taxass$abunds <- taxass$abunds[ ,index]
  return(taxass)
}

get.key.passing <- function(key, taxass){
  index <- is.element(el = key$limony.names, set = names(taxass$res))
  key <- key[index, ]
  return(key)
}

remove.empty.otus <- function(taxass){
  # summary(rowSums(taxass$abunds))
  index <- which(rowSums(taxass$abunds) == 0)
  cat("Removing", length(index), "OTUs that have zero reads.\n(this means they were only present in the controls or in samples that the total reads cutoff removed)\n")
  taxass$abunds <- taxass$abunds[-index, ]
  taxass$names <- taxass$names[-index, ]
  taxass$pvalues <- taxass$pvalues[-index, ]
  return(taxass)
}

add.pf.to.key.SampleName <- function(key, last.pf.date = "08Sep2003"){
  last.pf.date <- parse_date_time(x = last.pf.date, orders = "dmy")
  index <- order(key$Sample.Dates)
  key <- key[index, ]
  
  index.old <- which(key$Sample.Dates <= last.pf.date)
  index.ww <- which(key$Biological.Replicate == "W") # they all are distinguished with a W that I added manually in data entry, all uppercase
  index.pf <- setdiff(x = index.old, y = index.ww)
  
  key$Sample.Name[index.pf] <- paste0(key$Sample.Name[index.pf], "pf")
  
  # check
  # key[index.old, c(1:5,31)]
  return(key)
}

make.NA.reps.be.text <- function(key, column){
  index.c <- which(colnames(key) == column)
  index.r <- which(is.na(key[ ,index.c]))
  key[index.r, index.c] <- "NA"
  return(key)
}

get.change.from.info <- function(change.from, key, combinations){
  cat("\n\n\nDoing biol rep =",change.from, "\n")
  cat("All the filter rep options for this biol rep name: \n\n")
  index <- which(combinations$Biological.Replicate == change.from)
  print(combinations[index, ])
  
  row.indices <- which(key$Biological.Replicate == change.from)
  return(row.indices)
}

change.rep.to <- function(row.indices, change.to, key = key){
  those.samples <- key[row.indices, 1, drop = F]
  that.day <- merge(x = those.samples, y = key, by = "Sample.Name", all.x = T, all.y = F)
  day.combinations <- unique(that.day[ ,3:4])
  
  cat("\n\nAll the biol rep + filter rep combinations on days with this biol rep:\n\n")
  print(day.combinations)
  
  key$Bio.Rep[row.indices] <- change.to
  cat("\nChanged", length(row.indices), key$Biological.Replicate[row.indices[1]] ,"to", change.to, "\n")
  return(key)
}

change.biol.reps.to.letters <- function(key){
  key <- make.NA.reps.be.text(key = key, column = "Biological.Replicate")
  
  combinations <- unique(key[ ,c(3:4)])
  index <- order(combinations$Biological.Replicate, combinations$Filter.Replicate)
  combinations <- combinations[index, ]
  
  key <- cbind(key, "Bio.Rep" = key$Biological.Replicate)
  
  # All unknown replicates called the same biol replicate:
  
  index <- get.change.from.info(change.from = "NA", key = key, combinations = combinations)
  key <- change.rep.to(row.indices = index, change.to = "A", key = key)
  
  index <- get.change.from.info(change.from = "not recorded", key = key, combinations = combinations)
  key <- change.rep.to(row.indices = index, change.to = "A", key = key)
  
  index <- get.change.from.info(change.from = "unknown", key = key, combinations = combinations)
  key <- change.rep.to(row.indices = index, change.to = "A", key = key)
  
  # All bottle 1, Rep 1, and 1 to A
  
  index <- get.change.from.info(change.from = "B1", key = key, combinations = combinations)
  key <- change.rep.to(row.indices = index, change.to = "A", key = key)
  
  index <- get.change.from.info(change.from = "Rep1", key = key, combinations = combinations)
  key <- change.rep.to(row.indices = index, change.to = "A", key = key)
  
  index <- get.change.from.info(change.from = "Rep 1", key = key, combinations = combinations)
  key <- change.rep.to(row.indices = index, change.to = "A", key = key)
  
  index <- get.change.from.info(change.from = 1, key = key, combinations = combinations)
  key <- change.rep.to(row.indices = index, change.to = "A", key = key)
  
  # All W to A
  
  index <- get.change.from.info(change.from = "W", key = key, combinations = combinations)
  key <- change.rep.to(row.indices = index, change.to = "A", key = key)
  
  # All bottle 2, rep 2, and 2 to B
  
  index <- get.change.from.info(change.from = "B2", key = key, combinations = combinations)
  key <- change.rep.to(row.indices = index, change.to = "B", key = key)
  
  index <- get.change.from.info(change.from = "Rep2", key = key, combinations = combinations)
  key <- change.rep.to(row.indices = index, change.to = "B", key = key)
  
  index <- get.change.from.info(change.from = "Rep 2", key = key, combinations = combinations)
  key <- change.rep.to(row.indices = index, change.to = "B", key = key)
  
  index <- get.change.from.info(change.from = 2, key = key, combinations = combinations)
  key <- change.rep.to(row.indices = index, change.to = "B", key = key)
  
  # Change all B3 to C amd B4 to D
  index <- get.change.from.info(change.from = "B3", key = key, combinations = combinations)
  key <- change.rep.to(row.indices = index, change.to = "C", key = key) # The 1 rep instead of B1 rep is prob b/c I subbed in a different tube bc original tube was half- not a problem here
  
  index <- get.change.from.info(change.from = "B4", key = key, combinations = combinations)
  key <- change.rep.to(row.indices = index, change.to = "D", key = key)
  
  combinations <- unique(key[ ,c(37,4)])
  index <- order(combinations$Bio.Rep, combinations$Filter.Replicate)
  combinations <- combinations[index, ]
  
  return(key)
}

compress.letters <- function(key){
  all.samplenames <- unique(key$Sample.Name) 
  for(s in all.samplenames){
    index <- which(key$Sample.Name == s)
    days.reps <- key$Bio.Rep[index]
    days.reps <- unique(days.reps)
    days.reps <- sort(days.reps)
    
    if(days.reps[1] != "A"){
      index <- which(key$Sample.Name == s & key$Bio.Rep == days.reps[1])
      key$Bio.Rep[index] <- "A"
      cat("replaced", days.reps[1], "with A because it was the 1st rep on",s,"\n")
    }
    
    if (length(days.reps) > 1){
      if(days.reps[2] != "B"){
        index <- which(key$Sample.Name == s & key$Bio.Rep == days.reps[2])
        key$Bio.Rep[index] <- "B"
        cat("replaced", days.reps[2], "with B because it was the 2nd rep on",s,"\n")
      }
    }
    
    if (length(days.reps) > 2){
      if(days.reps[3] != "C"){
        index <- which(key$Sample.Name == s & key$Bio.Rep == days.reps[3])
        key$Bio.Rep[index] <- "C"
        cat("replaced", days.reps[3], "with C because it was the 3rd rep on",s,"\n")
      }
    }
    
    if (length(days.reps) > 2){
      cat("\n\n\n\nWOAH THIS EXISTS\n\n\n\n")
    }
  }
  return(key)
}

standardize.filter.numbers <- function(key){
  # Each unique filter code corresponds to a filter rep for that sample.
  # can just ignore the recorded "filter replicate" and use the filter codes
  # The filter codes also follow the same sequential order
  
  key <- make.NA.reps.be.text(key = key, column = "Filter.Replicate")
  key <- cbind(key, "Filt.Rep" = key$Filter.Replicate)
  
  all.samplenames <- unique(key$Sample.Name) 
  for(s in all.samplenames){
    index <- which(key$Sample.Name == s)
    filters <- key$Filter.Code[index]
    filters <- unique(filters)
    filters <- sort(filters)
    
    for (f in 1:length(filters)){
      index <- which(key$Filter.Code == filters[f])
      key$Filt.Rep[index] <- as.character(f)
    }
  }
  return(key)
}

add.extraction.replicates <- function(key){
  # Each unique extraction code corresponds to an extraction rep for that filter
  # for example, two halves of a filter may be extracted separately as rr 25 and rd 25
  # but I think I didn't sequence them both in that case, so that's why these are all rep 1 here
  
  key <- cbind(key, "Ext.Rep" = "")
  
  all.filternames <- unique(key$Filter.Code) 
  for(f in all.filternames){
    index <- which(key$Filter.Code == f)
    extractions <- key$Extraction.Code[index]
    extractions <- unique(extractions)
    extractions <- sort(extractions)
    
    for (e in 1:length(extractions)){
      index <- which(key$Extraction.Code == extractions[e])
      key$Ext.Rep[index] <- as.character(e)
    }
  }
  return(key)
}

add.sequencing.replicates <- function(key){
  # Each unique limony sequencing code corresponds to a sequencing rep for that extraction
  # for example, if there were low reads per sample I would try to re-sequence it on the next plate.
  
  key <- cbind(key, "Seq.Rep" = "")
  
  all.extractions <- unique(key$Extraction.Code) 
  for (e in all.extractions){
    index <- which(key$Extraction.Code == e)
    seqs <- key$limony.names[index]
    seqs <- unique(seqs)
    plate <- sub(pattern = "^.*_p", replacement = "", x = seqs)
    index <- order(as.numeric(plate))
    seqs <- seqs[index]
    
    for (s in 1:length(seqs)){
      index <- which(key$limony.names == seqs[s])
      key$Seq.Rep[index] <- as.character(s)
    }
  }
  return(key)
}

get.new.sample.colnames <- function(key){
  no.MEs <- sub(pattern = "^ME", replacement = "", x = key$Sample.Name)
  key <- cbind(key, "in.R.colnames" = no.MEs)
  
  # make different samples different biological reps once combined
  existing.Bio.Reps <- unique(key$Bio.Rep)
  existing.Bio.Reps <- sort(existing.Bio.Reps)
  highest.Bio.Rep <- which(LETTERS == existing.Bio.Reps[length(existing.Bio.Reps)])
  key <- cbind(key, "temp.bio.rep" = key$Bio.Rep)
  
  ext <- substr(x = key$Sample.Name, start = 12, stop = 100)
  index.pf <- grep(pattern = "pf", x = ext, value = F)
  ext <- sub(pattern = "pf", replacement = "", x = ext)
  index <- which(ext != "")
  ext.dates <- key$Sample.Dates[index]
  ext.dates <- unique(ext.dates)
  
  # add a dot before the pf in sample names with pf that match and won't change otherwise
  index.pf <- setdiff(x = index.pf, y = index)
  key$in.R.colnames[index.pf] <- sub(pattern = "pf", replacement = ".pf", x = key$in.R.colnames[index.pf])
  
  for (d in ext.dates){
    index <- which(key$Sample.Dates == d)
    pf.index <- grepl(pattern = "pf", x = key$Sample.Name[index])
    index.reg <- index[!pf.index]
    index.pf <- index[pf.index]
    indexes <- list("reg" = index.reg, "pf" = index.pf)
    
    for(i in 1:length(indexes)){
      if (length(indexes[[i]]) > 0){
        old.names <- key$Sample.Name[indexes[[i]]]
        new.name.base <- substr(x = old.names, start = 3, stop = 11)
        new.name.suffix <- substr(x = old.names, start = 12, stop = 100)
        new.name.suffix <- sub(pattern = "pf", replacement = "", x = new.name.suffix)
        new.name.suffix <- sub(pattern = "S", replacement = "D0", x = new.name.suffix)
        new.name.depth.suffix <- sub(pattern = "^.*D", replacement = "D", x = new.name.suffix)
        new.name.space.suffix <- sub(pattern = "D.*$", replacement = "", x = new.name.suffix)
        
        if (length(unique(new.name.base)) > 1){
          cat("\n\n   Shit why are different base sample names being combined???? \n\n")
          cat("on date",d,"combining",old.names,"\n\n")
        }
        
        if (length(unique(old.names)) > 1){ # specify that biol reps of diff sample names are different
          diff.samples <- unique(old.names)
          for (o in 1:length(diff.samples)){
            index <- which(key$Sample.Name == old.names[o])
            key$temp.bio.rep[index] <- LETTERS[highest.Bio.Rep + o]
            cat("Changed temp bio rep of", old.names[o], "from", key$Bio.Rep[index], "to", key$temp.bio.rep[index], "\n")
          }
        }
        
        if (length(unique(new.name.depth.suffix)) > 1){
          cat("Combining samples of different depth!!\n")
          new.name.depth.suffix <- "D"
        }
        
        if (length(unique(new.name.space.suffix)) > 1){
          cat("Combining samples with different spatial locations!!\n")
          new.name.space.suffix <- "s"
        }
        
        new.name <- paste0(new.name.base[1], ".", new.name.space.suffix[1], new.name.depth.suffix[1])
        if (names(indexes)[i] == "pf"){
          new.name <- paste0(new.name,"pf")
        }
        
        old.names <- sub(pattern = "^ME", replacement = "", x = old.names)
        key$in.R.colnames[indexes[[i]]] <- new.name
        cat(old.names, "    -->    ", new.name, "\n")
      }
    }
  }
  return(key)
}

assign.Lim.Rep <- function(key){
  lim.rep.key <- data.frame("temp.bio.rep" = key$temp.bio.rep, "Filt.Rep" = key$Filt.Rep, "Ext.Rep" = key$Ext.Rep)
  lim.rep.key <- unique(lim.rep.key)
  index <- order(lim.rep.key$Ext.Rep, lim.rep.key$Filt.Rep, lim.rep.key$temp.bio.rep)
  lim.rep.key <- lim.rep.key[index, ]
  lim.rep.key <- data.frame("temp" = paste(lim.rep.key$temp.bio.rep, lim.rep.key$Filt.Rep, lim.rep.key$Ext.Rep),   "Lim.Rep" = LETTERS[1:nrow(lim.rep.key)])
  
  key <- cbind(key, "temp" = paste(key$temp.bio.rep, key$Filt.Rep, key$Ext.Rep))
  
  key <- merge(x = lim.rep.key, y = key, by = "temp")
  unique(key$Lim.Rep)
  
  sample.names <- unique(key$in.R.colnames)
  for (s in sample.names){
    index.key <- which(key$in.R.colnames == s)
    reps <- key$Lim.Rep[index.key]
    reps <- unique(reps)
    reps <- sort(reps)
    
    for (r in 1:length(reps)){
      index.reps <- which(key$Lim.Rep[index.key] == reps[r])
      key$temp[index.key][index.reps] <- LETTERS[r]
    }
  }
  
  key$Lim.Rep <- key$temp
  index <- which(colnames(key) == "temp")
  key <- key[ ,-index]
  
  index <- which(colnames(key) == "temp.bio.rep")
  key <- key[ ,-index]
  
  index <- which(colnames(key) == "in.R.colnames")
  key <- key[ ,c(index, (1:ncol(key))[-index])]
  
  index <- order(key$Sample.Dates, key$Lim.Rep)
  key <- key[index, ]
  
  return(key)
}

# ---- pull out all good lake.samples ----

data.entry.key <- readRDS(file = file.path(extractions.data.entry.folder, input.data.entry.extractions))
taxass.list <- readRDS(file = file.path(current.r.data.folder, input.formatted.taxass))

key.samples <- get.key.samples(key = data.entry.key)
taxass.samples <- get.taxass.samples(taxass = taxass.list, key = key.samples)

taxass.samples <- get.taxass.passing(taxass = taxass.samples, cutoff = CHOSEN.CUTOFF)
key.samples <- get.key.passing(key = key.samples, taxass = taxass.samples)

taxass.samples <- remove.empty.otus(taxass = taxass.samples)


# ---- add pf to Sample.Name's that were prefiltered ----

key.samples <- add.pf.to.key.SampleName(key = key.samples)


# ---- assign consistent replicate names in key ----

# Biological Replicates are different Sample.Names (now that pf was added, counting WW as different)
  # named A, B, C, D ...
# Filter Replicates are different Filter.Codes (base on filter.code not filtering.replicate)
  # named 1, 2, 3, 4 ...
# Extraction Replicates are different Extraction.Codes (based on the rr99 code)  
  # named 1, 2, 3, 4 ...
# Sequencing Replicates are different limony.names (based on the rr99_p1 code)
  # named 1, 2, 3, 4 ...

key.samples <- change.biol.reps.to.letters(key = key.samples)

key.samples <- compress.letters(key = key.samples)

key.samples <- standardize.filter.numbers(key = key.samples)

key.samples <- add.extraction.replicates(key = key.samples)

key.samples <- add.sequencing.replicates(key = key.samples)

unique(key.samples[ ,3:4])
unique(key.samples[order(key.samples$Bio.Rep, key.samples$Filt.Rep, key.samples$Ext.Rep, key.samples$Seq.Rep),37:40])


# ---- Choose "limony" replicate ----

# For the time series analysis, everything gets grouped into either a biological replicate or a combined sample.
  # some diff names should be the same sample for big picture time-series analysis purposes. ex: 
    # winter S1D0 and S2D0 nearby stations
    # 2007 D7 and D12 two depths on same day
    # however WW and pf on same day should stay separate- may want to remove pf in some analyses (Expect bigger impact from pf than depth/station abnormal)
  # some technical reps will be considered "biological" reps in the limony dataset. ex: 
    # filter replicates ==> limony replicates (2 filters from same water grab)
    # extraction reps ==> limony replicates (2 extractions from diff halves of the filter) **THIS NEVER HAPPENS ANYWAY THOUGH**
  # some technical reps will be "pooled" for the limony analysis  
    # sequencing reps ==> pool into 1 sample with deeper sequencing depth (same DNA sequenced twice) (note both runs must pass reads/sample threshold)

key.samples <- get.new.sample.colnames(key = key.samples)

key.samples <- assign.Lim.Rep(key = key.samples)


# ---- export ----

my.file <- file.path(current.r.data.folder, created.samples.key)
cat("Making file: ", my.file, "\n")
# saveRDS(object = key.samples, file = my.file)

my.file <- file.path(current.r.data.folder, created.samples.list)
cat("Making file: ", my.file, "\n")
# saveRDS(object = taxass.samples, file = my.file)

