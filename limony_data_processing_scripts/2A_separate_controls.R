# RRR 7/12/20
# remove the controls (b/c many blanks "fail" sequencing)
# These will be dealt with separately from the main analyses anyway
# So save them in their own list and with their own key


# ---- set-up ----

input.formatted.taxass <- "data/limony-IRD/2021-08-25_processing/1B_taxass.rds"
input.data.entry.extractions <- "data/metadata/2020-07-10_data_entry_limony_and_TYMEFLIES.rds"

output.folder <- "data/limony-IRD/2021-08-25_processing"

created.controls.list <- "2A_taxass_controls.rds"
created.controls.key <- "2A_key_controls.rds"

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

get.key.controls <- function(key){
  index <- which(key$sample.type != "lake.sample" & !is.na(key$limony.names))
  key <- key[index, ]
  
  kit.number <- sub(pattern = "kit.", replacement = "", x = key$kit.number)
  kit.number <- as.numeric(kit.number)
  index <- order(key$sample.type, key$Sample.Dates, key$Filter.Code, kit.number, key$Filter.Replicate)
  key <- key[index, ]
  rownames(key) <- 1:nrow(key)
  
  return(key)
}

get.taxass.controls <- function(taxass, key){
  limony.names <- names(taxass$res)
  
  index <- is.element(el = limony.names, set = key$limony.names)
  taxass$res <- taxass$res[index]
  taxass$abunds <- taxass$abunds[ ,index]
  
  cat(all.equal(names(taxass$res), colnames(taxass$abunds)))
  return(taxass)
}

create.colnames.for.controls <- function(key){
  index.blanks <- which(key$sample.type == "blank")
  blanks <- key[index.blanks, ]
  blank.nums <- sub(pattern = "kit\\.\\d*\\.", replacement = "", x = blanks$Filter.Code) 
  blanks.kit <- sub(pattern = "\\.blank.*$", replacement = "", x = blanks$Filter.Code)
  blanks.plate <- sub(pattern = "^plate\\.", replacement = "", x = blanks$limony.plate)
  blanks.plate <- sub(pattern = "\\..*$", replacement = "", x = blanks.plate)
  blanks <- cbind("in.r.colnames" = paste0(blanks.kit, ".", blank.nums, ".p", blanks.plate),
                  blanks)
  
  index.gendonors <- which(key$sample.type == "generous.donor")
  gendonors <- key[index.gendonors, ]
  gen.date <- sub(pattern = "ME", replacement = "", x = gendonors$Sample.Name)
  gen.plate <- sub(pattern = "^plate\\.", replacement = "", x = gendonors$limony.plate)
  gen.plate <- sub(pattern = "\\..*$", replacement = "", x = gen.plate)
  gendonors <- cbind("in.r.colnames" = paste0("gendonor.", gen.date, ".p", gen.plate),
                     gendonors)
  
  index.mocks <- which(key$sample.type == "mock")
  mocks <- key[index.mocks, ]
  mock.type <- sub(pattern = "^.*mock.", replacement = "", x = mocks$Filter.Code)
  mock.plate <- sub(pattern = "^plate\\.", replacement = "", x = mocks$limony.plate)
  mock.plate <- sub(pattern = "\\..*$", replacement = "", x = mock.plate)
  mocks <- cbind("in.r.colnames" = paste0("mock.", mock.type, ".p", mock.plate),
                 mocks)
  
  key <- rbind(blanks, gendonors, mocks)
  return(key)
}

name.taxass.controls <- function(taxass, key){
  taxass.names <- data.frame("limony.names" = colnames(taxass$abunds), "taxass.order" = 1:ncol(taxass$abunds))
  key <- cbind(key, "key.order" = as.numeric(row.names(key)))
  key <- merge(x = taxass.names, y = key, by = "limony.names", all = F) # watch, nrow stays the same good
  
  # everything in order of the taxass objects
  index <- order(key$taxass.order)
  key <- key[index, ]
  cat(all.equal(colnames(taxass$abunds), names(taxass$res)))
  cat(all.equal(colnames(taxass$abunds), key$limony.names))
  
  # everything in order of the organized key
  index <- order(key$key.order)
  taxass$abunds <- taxass$abunds[ ,index]
  taxass$res <- taxass$res[index]
  key <- key[index, ]
  cat(all.equal(colnames(taxass$abunds), names(taxass$res)))
  cat(all.equal(colnames(taxass$abunds), key$limony.names))
  
  # rename taxass colnames from the matching key
  colnames(taxass$abunds) <- key$in.r.colnames
  names(taxass$res) <- key$in.r.colnames
  cat(all.equal(colnames(taxass$abunds), names(taxass$res)))
  
  return(taxass)
}

# ---- pull out controls ----

data.entry.key <- readRDS(file = input.data.entry.extractions)
key.ctrls <- get.key.controls(key = data.entry.key)

taxass.list <- readRDS(file = input.formatted.taxass)
taxass.ctrls <- get.taxass.controls(taxass = taxass.list, key = key.ctrls)

key.ctrls <- create.colnames.for.controls(key = key.ctrls)
taxass.ctrls <- name.taxass.controls(taxass = taxass.ctrls, key = key.ctrls)

# ---- export ----

my.file <- file.path(output.folder, created.controls.list)
cat("Making file: ", my.file, "\n")
# saveRDS(object = taxass.ctrls, file = my.file)

my.file <- file.path(output.folder, created.controls.key)
cat("Making file: ", my.file, "\n")
# saveRDS(object = key.ctrls, file = my.file)
