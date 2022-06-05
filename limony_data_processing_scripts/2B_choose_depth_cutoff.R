# RRR 7/12/20
# Choose the sequencing depth cutoff for what should be included
# and what should be considered a failed run and discarded

# ---- set-up ----

input.formatted.taxass <- "data/limony-IRD/2021-08-25_processing/1B_taxass.rds"
input.data.entry.extractions <- "data/metadata/2020-07-10_data_entry_limony_and_TYMEFLIES.rds"

plots.folder <- "plots/2021-08-25_Limony_IRD_seq_depth_cutoff_choice"

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

get.samples.only <- function(key, taxass){
  index <- which(key$sample.type == "lake.sample" & !is.na(key$limony.names))
  key <- key[index, ]
  
  index <- is.element(el = names(taxass$res), set = key$limony.names)
  taxass$res <- taxass$res[index]
  taxass$abunds <- taxass$abunds[ ,index]
  
  cat(all.equal(names(taxass$res), colnames(taxass$abunds)))
  return(taxass)
}

plot.reads.per.sample.line <- function(res = taxass.samples$res, cutoff = "none", x.lim = c(1, length(res)), y.lim = c(min(res), max(res))){
  res <- sort(res)
  res.stats <- summary(res)
  
  point.cols <- rep.int(x = adjustcolor(col = "black", alpha.f = .3), times = length(res))
  if (cutoff != "none"){
    index <- res < cutoff
    point.cols[index] <- adjustcolor(col = "orange", alpha.f = .5)
    num.removed <- sum(index)
  }
  
  par(mar = c(3.5,6.5,.5,.5))
  plot(x = x.lim, y = y.lim, type = "n", xlim = x.lim, ylim = y.lim, axes = F, ann = F)
  points(x = 1:length(res), y = res, pch = 16, col = point.cols, bg = point.cols, cex = .5)
  axis(side = 1, labels = F, at = x.lim)
  axis(side = 1, tick = F, at = x.lim, line = -.5)
  axis(side = 2, at = c(y.lim,res.stats[2], res.stats[3], res.stats[5]), labels = F)
  axis(side = 2, at = c(y.lim,res.stats[2], res.stats[3], res.stats[5]), labels = c(y.lim,"Q1","Med","Q3"), tick = F, las = 1, line = -.25)
  # mtext(text = "Reads per Sample", side = 3, line = .5)
  mtext(text = "Sample (rank order)", side = 1, line = 1)
  mtext(text = "Reads", side = 2, line = 4, at = y.lim[2] - y.lim[2] / 5)
  
  if (cutoff != "none"){
    # lines(x = c(-50, num.removed), y = c(cutoff, cutoff), col = "orange")
    # lines(x = c(num.removed, num.removed), y = c(-500, cutoff), col = "orange", xpd = T)
    axis(side = 1, at = num.removed, col.axis = "orange", col.ticks = "orange", tck = -.1, lwd = 3, label = F)
    axis(side = 1, at = num.removed, col.axis = "orange", col.ticks = "orange", tick = F, label = T, line = 1.3, cex.axis = 1.5)
    axis(side = 2, at = cutoff, col.axis = "orange", col.ticks = "orange", tck = -.12, lwd = 3, label = F)
    axis(side = 2, at = cutoff, col.axis = "orange", col.ticks = "orange", tick = F, label = T, line = 1.5, las = 2, cex.axis = 1.5)
    
  }
  
}

plot.reads.per.sample.hist <- function(res = taxass.samples$res, cutoff = "none"){
  x.lim = c(1, length(res))
  y.lim = c(min(res), max(res))
  
  h <- hist(x = res, breaks = 100, plot = F)
  d <- density(x = res)
  x.ax <- seq(from = min(res), to = max(res), length.out = 1000)
  n <- dnorm(x = x.ax, mean = mean(res), sd = sd(res))
  
  bin.cols <- rep.int(x = adjustcolor(col = "black", alpha.f = .3), times = length(h$breaks))
  if (cutoff != "none"){
    index <- h$breaks < cutoff
    bin.cols[index] <- adjustcolor(col = "orange", alpha.f = .5)
  } 
  
  par(mar = c(4,4,2,.5))
  plot(h, freq = F, ann = F, col = bin.cols)
  lines(d, lwd = 3, col = adjustcolor(col = "black", alpha.f = .5))
  lines(x = x.ax, y = n, col = adjustcolor(col = "navy", alpha.f = .3), lwd = 1)
  # mtext(text = "Histogram of Sequencing Depth", side = 3, line = .5)
  mtext(text = "Reads per Sample", side = 1, line = 2.5)
  mtext(text = "Density", side = 2, line = 2.5)
  mtext(text = "- normal curve", side = 1, line = -4, at = 50000, col = adjustcolor(col = "navy", alpha.f = .3), cex = .7, adj = 0)
  mtext(text = "- sample density", side = 1, line = -3, at = 50000, col = adjustcolor(col = "black", alpha.f = .5), cex = .7, adj = 0)
  
  if (cutoff != "none"){
    abline(v = cutoff, col = "orange")
    mtext(text = cutoff, side = 3, at = cutoff, adj = 1, col = "orange")
  }
}  

plot.reads.per.sample.cloud <- function(res = taxass.samples$res, cutoff = "none", x.lim = c(0, max(res))){
  res <- sort(res)
  
  point.cols <- rep.int(x = adjustcolor(col = "black", alpha.f = .3), times = length(res))
  # if (cutoff != "none"){
  #   index <- res < cutoff
  #   point.cols[index] <- adjustcolor(col = "orange", alpha.f = .5)
  #   num.removed <- sum(index)
  # }
  
  y.vals <- runif(n = length(res), min = 0, max = 1)
  
  par(mar = c(4,2,2,.5))
  plot(x = x.lim, y = c(0,1), type = "n", xlim = x.lim, ylim = c(0,1), ann = F, axes = F)
  points(x = res, y = y.vals, pch = 16, col = point.cols, bg = point.cols)
  box(which = "plot", col = adjustcolor("black",.5))
  axis(side = 1)
  mtext(text = "Reads per Sample", side = 1, line = 2.5)
  mtext(text = "each dot is a sample", side = 2, line = .5)
  
  if (cutoff != "none"){
    abline(v = cutoff, col = "orange")
    mtext(text = cutoff, side = 3, at = cutoff, adj = 1, col = "orange")
  }
  
}


# ---- analysis ----

taxass.list <- readRDS(file = input.formatted.taxass)
data.entry.key <- readRDS(file = input.data.entry.extractions)

taxass.samples <- get.samples.only(key = data.entry.key, taxass = taxass.list)

# ---- make plots ----

# no cutoff
plot.reads.per.sample.line(res = taxass.samples$res)
plot.reads.per.sample.line(res = taxass.samples$res, y.lim = c(0,20000), x.lim = c(1,200))
plot.reads.per.sample.hist(res = taxass.samples$res)
plot.reads.per.sample.cloud(res = taxass.samples$res)
plot.reads.per.sample.cloud(res = taxass.samples$res, x.lim = c(0,30000))

# yes cutoff
c = 8500
plot.reads.per.sample.line(res = taxass.samples$res, cutoff = c)
plot.reads.per.sample.line(res = taxass.samples$res, y.lim = c(0,30000),  cutoff = c)
plot.reads.per.sample.line(res = taxass.samples$res, y.lim = c(0,20000), x.lim = c(1,500), cutoff = c)
plot.reads.per.sample.hist(res = taxass.samples$res, cutoff = c)
plot.reads.per.sample.cloud(res = taxass.samples$res, cutoff = c)
plot.reads.per.sample.cloud(res = taxass.samples$res, x.lim = c(0,30000), cutoff = c)


m = matrix(nrow = 3, ncol = 2, data = c(1,1,2,4,3,4), byrow = T)
layout(mat = m, widths = c(1.5,1), heights = c(1.5,1,1))
par(oma = c(.25,.25,3,.25))
plot.reads.per.sample.line(res = taxass.samples$res, y.lim = c(0,20000), x.lim = c(1,500), cutoff = c)
mtext(text = "Sequencing Depth of limony Samples", side = 3, line = 1, outer = T, cex = 1.5)
plot.reads.per.sample.line(res = taxass.samples$res, cutoff = c)
plot.reads.per.sample.hist(res = taxass.samples$res, cutoff = c)
plot.reads.per.sample.cloud(res = taxass.samples$res, x.lim = c(0,30000), cutoff = c)

box(which = "inner", col="red", lwd = 3)
box(which = "outer", col="blue", lwd = 3)
box(which = "plot", col="purple", lwd = 3)
box(which = "figure", col="orange", lwd = 3)

# ---- export plots ----

my.plot <- file.path(plots.folder, "no_cutoff.pdf")
cat("Making file: ", my.plot, "\n")
pdf(file = my.plot, width = 10, height = 6.5)

m = matrix(nrow = 3, ncol = 2, data = c(1,1,2,4,3,4), byrow = T)
layout(mat = m, widths = c(1.5,1), heights = c(1.5,1,1))
par(oma = c(.25,.25,3,.25))
plot.reads.per.sample.line(res = taxass.samples$res, y.lim = c(0,20000), x.lim = c(1,500))
mtext(text = "Sequencing Depth of limony Samples", side = 3, line = 1, outer = T, cex = 1.5)
plot.reads.per.sample.line(res = taxass.samples$res)
plot.reads.per.sample.hist(res = taxass.samples$res)
plot.reads.per.sample.cloud(res = taxass.samples$res, x.lim = c(0,30000))

dev.off()


c = 10000
my.plot <- file.path(plots.folder, "cutoff = 10,000.pdf")
cat("Making file: ", my.plot, "\n")
pdf(file = my.plot, width = 10, height = 6.5)

m = matrix(nrow = 3, ncol = 2, data = c(1,1,2,4,3,4), byrow = T)
layout(mat = m, widths = c(1.5,1), heights = c(1.5,1,1))
par(oma = c(.25,.25,3,.25))
plot.reads.per.sample.line(res = taxass.samples$res, y.lim = c(0,20000), x.lim = c(1,500), cutoff = c)
mtext(text = "Sequencing Depth of limony Samples", side = 3, line = 1, outer = T, cex = 1.5)
plot.reads.per.sample.line(res = taxass.samples$res, cutoff = c)
plot.reads.per.sample.hist(res = taxass.samples$res, cutoff = c)
plot.reads.per.sample.cloud(res = taxass.samples$res, x.lim = c(0,30000), cutoff = c)
mtext(text = paste0("cutoff = 10,000\n1 read = ", round(x = 1/c * 100, digits = 3)," %"), side = 3, line = -15, outer = T, at = .7, col = "orange", adj = 0)

dev.off()



c = 9000
my.plot <- file.path(plots.folder, "cutoff = 9,000.pdf")
cat("Making file: ", my.plot, "\n")
pdf(file = my.plot, width = 10, height = 6.5)

m = matrix(nrow = 3, ncol = 2, data = c(1,1,2,4,3,4), byrow = T)
layout(mat = m, widths = c(1.5,1), heights = c(1.5,1,1))
par(oma = c(.25,.25,3,.25))
plot.reads.per.sample.line(res = taxass.samples$res, y.lim = c(0,20000), x.lim = c(1,500), cutoff = c)
mtext(text = "Sequencing Depth of limony Samples", side = 3, line = 1, outer = T, cex = 1.5)
plot.reads.per.sample.line(res = taxass.samples$res, cutoff = c)
plot.reads.per.sample.hist(res = taxass.samples$res, cutoff = c)
plot.reads.per.sample.cloud(res = taxass.samples$res, x.lim = c(0,30000), cutoff = c)
mtext(text = paste0("cutoff = 9,000\n1 read = ", round(x = 1/c * 100, digits = 3)," %"), side = 3, line = -15, outer = T, at = .7, col = "orange", adj = 0)

dev.off()



c = 8500
my.plot <- file.path(plots.folder, "cutoff = 8,500.pdf")
cat("Making file: ", my.plot, "\n")
pdf(file = my.plot, width = 10, height = 6.5)

m = matrix(nrow = 3, ncol = 2, data = c(1,1,2,4,3,4), byrow = T)
layout(mat = m, widths = c(1.5,1), heights = c(1.5,1,1))
par(oma = c(.25,.25,3,.25))
plot.reads.per.sample.line(res = taxass.samples$res, y.lim = c(0,20000), x.lim = c(1,500), cutoff = c)
mtext(text = "Sequencing Depth of limony Samples", side = 3, line = 1, outer = T, cex = 1.5)
plot.reads.per.sample.line(res = taxass.samples$res, cutoff = c)
plot.reads.per.sample.hist(res = taxass.samples$res, cutoff = c)
plot.reads.per.sample.cloud(res = taxass.samples$res, x.lim = c(0,30000), cutoff = c)
mtext(text = paste0("cutoff = 8,500\n1 read = ", round(x = 1/c * 100, digits = 3)," %"), side = 3, line = -15, outer = T, at = .7, col = "orange", adj = 0)

dev.off()



c = 8000
my.plot <- file.path(plots.folder, "cutoff = 8,000.pdf")
cat("Making file: ", my.plot, "\n")
pdf(file = my.plot, width = 10, height = 6.5)

m = matrix(nrow = 3, ncol = 2, data = c(1,1,2,4,3,4), byrow = T)
layout(mat = m, widths = c(1.5,1), heights = c(1.5,1,1))
par(oma = c(.25,.25,3,.25))
plot.reads.per.sample.line(res = taxass.samples$res, y.lim = c(0,20000), x.lim = c(1,500), cutoff = c)
mtext(text = "Sequencing Depth of limony Samples", side = 3, line = 1, outer = T, cex = 1.5)
plot.reads.per.sample.line(res = taxass.samples$res, cutoff = c)
plot.reads.per.sample.hist(res = taxass.samples$res, cutoff = c)
plot.reads.per.sample.cloud(res = taxass.samples$res, x.lim = c(0,30000), cutoff = c)
mtext(text = paste0("cutoff = 8,000\n1 read = ", round(x = 1/c * 100, digits = 3)," %"), side = 3, line = -15, outer = T, at = .7, col = "orange", adj = 0)

dev.off()



c = 6000
my.plot <- file.path(plots.folder, "cutoff = 6,000.pdf")
cat("Making file: ", my.plot, "\n")
pdf(file = my.plot, width = 10, height = 6.5)

m = matrix(nrow = 3, ncol = 2, data = c(1,1,2,4,3,4), byrow = T)
layout(mat = m, widths = c(1.5,1), heights = c(1.5,1,1))
par(oma = c(.25,.25,3,.25))
plot.reads.per.sample.line(res = taxass.samples$res, y.lim = c(0,20000), x.lim = c(1,500), cutoff = c)
mtext(text = "Sequencing Depth of limony Samples", side = 3, line = 1, outer = T, cex = 1.5)
plot.reads.per.sample.line(res = taxass.samples$res, cutoff = c)
plot.reads.per.sample.hist(res = taxass.samples$res, cutoff = c)
plot.reads.per.sample.cloud(res = taxass.samples$res, x.lim = c(0,30000), cutoff = c)
mtext(text = paste0("cutoff = 6,000\n1 read = ", round(x = 1/c * 100, digits = 3)," %"), side = 3, line = -15, outer = T, at = .7, col = "orange", adj = 0)

dev.off()


# ~end~