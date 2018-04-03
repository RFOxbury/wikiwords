#!/usr/local/bin/Rscript

# TEST:
#
# ./make_shiny_data.R ../wikipedia ./tmp
# OR for testing:
# args <- c("../wikipedia","tmp")

args <- commandArgs(trailingOnly=T)

args <- c("./data", "./shiny/wordrep")

inpath <- args[1]
outpath <- args[2]

HOMEPATH <- "/Users/bill/Projects/201801-wikiwords/wikiwords"
source(sprintf("%s/tda_functions.R", HOMEPATH))

# CONVENTION 
# look in 'inpath' for 3 files 'wdict.txt', 'wcounts.txt', 'w.vectors.txt'
# write to 'outpath' the following:
#
# * copy of 'ui.R', 'server.R', 'helpers.R'
# * subdirectory ./data containing files
#     data/w2v_data.Rds
#     data/w2v_cset.Rds
#     data/w2v_graph.Rds
#     data/w2v_layout.Rds

# check output path and create directory if it isn't already there:

cat("Checking for output directory...\n")

tmp <- ifelse(!dir.exists(file.path(args[2])), 
              dir.create(file.path(args[2])), 
              FALSE)
tmp <- ifelse(!dir.exists(file.path(sprintf("%s/data", args[2]))), 
              dir.create(file.path(sprintf("%s/data", args[2]))), 
              FALSE)

# read input data:

cat("Reading data...\n")

dictfile <- sprintf("%s/wdict.txt", args[1])
freqfile <- sprintf("%s/wcounts.txt", args[1])
vecfile <- sprintf("%s/wvectors.txt", args[1])

dat1 <- read.table(dictfile)
names(dat1) <- c("word", "class")
dat2 <- read.table(freqfile, row.names=NULL)
names(dat2) <- c("freq", "word")
dat2$freq <- as.numeric(dat2$freq)
dat3 <- read.table(vecfile)
names(dat3)[1] <- "word"
dat4 <- merge(dat1, dat2)
alldat <- merge(dat3, dat4)
attach(alldat)

cat("Fitting frequency distribution...\n")

# frequency distribution:
tab <- table(freq)
xx <- as.numeric(dimnames(tab)$freq)
yy <- as.numeric(tab)

imgfile <- sprintf("%s/word_frequency_distribution.pdf", args[2])
pdf(imgfile)
plot(yy ~ xx, log='xy', col='blue', cex=0.5, frame.plot=0, 
     xlab="Word count", ylab="Nr words")

# fit power law:
cond <- (xx < 2e03)
model <- lm(log(yy[cond]) ~ log(xx[cond]))
( coef <- model$coefficients )
f <- function(x) exp(coef[1] + coef[2]*log(x))
xxx <- exp(0:10)
yyy <- sapply(xxx, f)
points( yyy~xxx, type='l', col='red', lwd=2)
text(2e03, 1000, 
     labels=sprintf("Power law coefficients (%g, %g)", coef[1], coef[2]), col='red')
dev.off()

cat(sprintf("Plot written to %s/word_frequency_distribution.pdf\n", args[2]))

################################################################
# exponential binning of frequency using the power law model:

nbins <- 20                 # PARAMETER for mapper construction
coarseness <- 32            # PARAMETER target average cluster size

bins <- makebins(nbins, coef)

################################################################
# clustering 

# hierarchical clustering within bins:

cat(sprintf("Building clusters (using %d bins and average cluster size %d)...\n", 
            nbins, coarseness))
cluster.set <- makeclusters_centroids(bins, coarseness)
cat(sprintf("Found %d clusters\n", length(cluster.set)))

################################################################
# compute graph edges:

library(igraph)

cat("Computing mapper graph...\n")
edges <- edgelist(cluster.set)

mg <- graph_from_edgelist(edges, directed=FALSE)
cc <- clusters(mg)
cat(sprintf("Found %d graph edges, %d components, sizes:\n", 
            length(unique(as.numeric(edges))),
            cc$no))
table(cc$csize)

################################################################
# 2-dim layout using t-SNE:

library(Rtsne)

cat("Computing t-SNE layout...\n")
centroid.mat <- do.call(rbind, lapply(cluster.set, function(x) x$centroid))
lout <- Rtsne(centroid.mat, 
              theta=0.1,
              verbose=TRUE, 
              initial_dims=100,
              check_duplicates=FALSE, 
              max_iter=2000)
lout <- data.frame(lout$Y)
names(lout) <- c("X","Y")


################################################################
# write output:

cat("Writing output...\n")

path <- sprintf("%s/data", outpath)
dat <- dat2

save(dat, file=sprintf("%s/w2v_freq.Rds", path))
save(cluster.set, file=sprintf("%s/w2v_cset.Rds", path))
save(mg, file=sprintf("%s/w2v_graph.Rds", path))
save(lout, file=sprintf("%s/w2v_layout.Rds", path))

system(sprintf("cp ./shiny/*.R %s", args[2]))

cat("Done.\n")
