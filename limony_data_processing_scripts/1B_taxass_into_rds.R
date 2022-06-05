# RRR 5/26/20 ----
# At the end of TaxAss your "data" folder will contain:
#   - your starting QC'ed abundances file, like dada's seqtab_nochim.rds or mothur's .count_table
#   - your starting .fasta file
#   - created (in the first taxass formatting step) .abund file (tab-delim normalized OTU table)
#   - created (in the first taxass formatting step) .fasta file (unaligned fasta format)
#   - created (in the first taxass formatting step) .count file (to keep track of total reads when normalize OTU table)
#   - your final .tax file (the whole point of running TaxAss, this is what it actually makes)
# This script reads in the .tax, .abund, and .count files into R
# Eventually it may be incorporated into taxass if I make more export options available
# Next step would be dataset-dependent formatting, for example of the sample names

# ---- input ----

folder.path <- "data/limony-IRD/2021-08-25_processing/2021-08-25_taxass_output"

count.table.file <- file.path(folder.path, "IRD.count")
file.type.of.count <- "dada"
taxonomy.file <- file.path(folder.path, "IRD.98.80.80.taxonomy")
abund.file <- file.path(folder.path, "IRD.abund")

output.folder <- "data/limony-IRD/2021-08-25_processing"
output.filename <- "1B_taxass.rds"

# ---- functions ----

import.taxonomy <- function(my.file){
  tax <- read.csv(file = my.file, header = T, colClasses = "character")
  index <- order(tax$seqID)
  tax <- tax[index, ]
  tax <- as.matrix(tax)
  tax <- tax[ ,c(2:ncol(tax),1)] # move seqID to right side
  row.names(tax) <- NULL
  return(tax)
}

import.abundances <- function(my.file){
  abund <- read.table(file = my.file, header = T, colClasses = "character", sep = "\t")
  index <- order(abund$seqID)
  abund <- abund[index, ]
  seqIDs <- abund[ ,1]
  abund <- as.matrix(abund[ ,-1]) # move seqID to row name
  abund <- apply(X = abund, MARGIN = 2, FUN = as.numeric)
  row.names(abund) <- seqIDs
  index <- order(colnames(abund))
  abund <- abund[ ,index]
  return(abund)
}

import.raw.abundances <- function(my.file, FileType){
  if (FileType == "mothur"){
    count.table <- read.table(file = my.file, header = TRUE, sep = "\t", colClasses = "character")
    count.table <- count.table[ ,1:2]
    colnames(count.table) <- c("Sample.Name", "Total.Reads")
  }else if (FileType == "dada"){
    count.table <- read.table(file = my.file, header = T, sep = "\t", colClasses = "character")
  }else{
    return(cat("file.type.of.count input must be either \"dada\" or \"mothur\""))
  }
  count.vector <- count.table$Total.Reads
  count.vector <- as.numeric(count.vector)
  names(count.vector) <- count.table$Sample.Name
  index <- order(names(count.vector))
  count.vector <- count.vector[index]
  return(count.vector)
}

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

remove.parentheses <- function(x){
  # Remove parentheses and bootstrap % confidence from taxonomy names
  # x is a single name, use this function with apply
  fixed.name <- sub(pattern = '\\(.*\\)' , replacement = '', x = x)
  return(fixed.name)
}

extract.bootstrap.pvalues <- function(x){
  pvalue <- sub(pattern = "^.*\\(", replacement = "", x = x)
  pvalue <- sub(pattern = "\\)$", replacement = "", x = pvalue)
  pvalue <- as.numeric(pvalue)
  return(pvalue)
}


# ---- go ----

count.vector <- import.raw.abundances(my.file = count.table.file, FileType = file.type.of.count)
tax.matrix <- import.taxonomy(my.file = taxonomy.file)
abund.matrix <- import.abundances(my.file = abund.file)

pvalue.matrix <- apply(X = tax.matrix, MARGIN = 2, FUN = extract.bootstrap.pvalues)
# ignore error messages about NA's- that's bc it becomes NA if no parentheses (ie if "unclassified" or if "otu_1")
row.names(pvalue.matrix) <- tax.matrix[ ,8]
tax.matrix <- apply(X = tax.matrix, MARGIN = 2, FUN = remove.parentheses)

top(tax.matrix)
top(pvalue.matrix)
top(X = count.vector)
top(abund.matrix)

all.equal(row.names(abund.matrix), tax.matrix[ ,8])
all.equal(row.names(pvalue.matrix), tax.matrix[ ,8])
all.equal(colnames(abund.matrix), names(count.vector))

taxass <- list("names" = tax.matrix, "abunds" = abund.matrix, "res" = count.vector, "pvalues" = pvalue.matrix)
str(taxass)

my.file <- file.path(output.folder, output.filename)
cat("Making File: ", my.file, "\n")
# saveRDS(object = taxass, file = my.file)
