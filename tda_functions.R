library(shiny)
library(igraph)
library(wordcloud)

findnodes <- function(w){
  if(!(w %in% word)) return()
  id <- which(word==w)
  which( sapply(cluster.set, function(s){id %in% s$cluster}) )
}


makebins <- function(nbins, coef){
  
  epsilon <- 1/nbins/5
  left <- (0:(nbins-1))/nbins - epsilon
  right <- (1:nbins)/nbins + epsilon
  # aa <- exp(coef[1]) # don't need
  bb <- -coef[2]
  mn <- min(freq)^(1-bb)
  mx <- max(freq)^(1-bb)
  func <- function(x) x^(1/(1-bb))
  left <- func( mn + (mx - mn)*left )
  right <- func( mn + (mx - mn)*right )
  right[nbins] <- 1e09
  
  bins <- list()
  for(i in 1:nbins)
    bins[[i]] <- which(freq > left[i] & freq < right[i]) 
  # prune:
  empties <- c()
  for(i in 1:nbins) if(length(bins[[i]])==0) empties <- c(empties, i)
  for(e in rev(empties)){ bins[[e]] <- NULL; nbins <<- nbins-1 }
  # return 
  bins
}


makeclusters_centroids <- function(bins, coarseness){
  
  cluster.set <- list()
  ccount <- 0
  
  recurse.func <- function(hc, bin.nr, ccount){
    
    for(child in 1:2){
      set <- unique(as.numeric(labels(hc[[child]])))
      if(length(set) < coarseness){
        # update global variables:
        ccount <- ccount + 1
        centroid <- colSums( alldat[set, 2:201] )/length(set)
        cluster.set[[ccount]] <<- list(cluster=set, centroid=centroid, height=bin.nr)
        cat(sprintf("%d %d       \r", length(set), ccount))
      } else {
        ccount <- recurse.func( hc[[child]], bin.nr, ccount )
      }
    }
    return(ccount)
  }
  
  nbins <- length(bins)
  for(nr in 1:nbins){
    
    bin <- bins[[nr]]
    cat(sprintf("Bin %d                           \n", nr))
    bindist <- dist(alldat[bin, 2:201], method="euclidean") 
    hc <- as.dendrogram( hclust(bindist) )
    ccount <- recurse.func(hc, nr, ccount)
    cat("\r")
  }  
  
  # return deduped set:
  unique(cluster.set)
}


# CANDIDATE FOR speed up with h2o:
edgelist <- function(cluster.set){
  
  ccount <- length(cluster.set)
  edges <- matrix(nrow=0, ncol=2)
  for(i in 1:(ccount-1))
    for(j in (i+1):ccount){
      si <- cluster.set[[i]]$cluster
      sj <- cluster.set[[j]]$cluster
      wt <- length(si) + length(sj) - length(unique(c(si,sj)))
      if(wt > 0)
        edges <- rbind(edges, c(i,j))
    }
  mx <- max(edges)
  if(mx < ccount)
    for(i in (mx+1):ccount)
      edges <- rbind(edges, c(i,i))
  # return:
  edges
}


components <- function(g){
  
  cc <- clusters(g)
  # return:
  lapply(1:cc$no, function(i){ 
    induced.subgraph(g, which(cc$membership == i)) })
}

compsets <- function(g, cset, j){
  
  cc <- clusters(g)
  set <- as.numeric( V(g)[cc$membership==j] )
  # return:
  lapply(set, function(s) cset[[s]])
}

# hand-craft the plot function:
handplot <- function(g, lout){
  
  mai <- par()$mai
  par(mai=c(0,0,0,0))
  plot(lout, type='n', axes=FALSE, xlab="", ylab="")
  for(e in 1:ecount(g)){
    v <- ends(g,e)
    segments(lout[v[1],1], lout[v[1],2],  lout[v[2],1], lout[v[2],2], col='grey')
  }
  points(lout, pch=1, cex=1.5*V(g)$size, col= "grey")
  points(lout, pch=19, cex=V(g)$size, col= V(g)$color)
  xtext <- min(lout[,1]) + 0.1*(max(lout[,1]) - min(lout[,1]))
  ytext <- min(lout[,2]) + 0.01*(max(lout[,2]) - min(lout[,2]))
  text(xtext, ytext, labels=sprintf("%d vertices, %d edges", vcount(g), ecount(g)), 
       col='darkblue')
  par(mai = mai)
}

handsubplot <- function(g, lout){
  
  mai <- par()$mai
  par(mai=c(0,0,0,0))
  plot(lout, type='n', axes=FALSE, xlab="", ylab="")
  for(e in 1:ecount(g)){
    v <- ends(g,e)
    segments(lout[v[1],1], lout[v[1],2],  lout[v[2],1], lout[v[2],2], col='grey')
  }
  points(lout, pch=1, cex=3.5*V(g)$size, col= "grey")
  points(lout, pch=19, cex=3*V(g)$size, col= V(g)$color)
  xtext <- min(lout[,1]) + 0.1*(max(lout[,1]) - min(lout[,1]))
  ytext <- min(lout[,2]) + 0.01*(max(lout[,2]) - min(lout[,2]))
  text(xtext, ytext, labels=sprintf("%d vertices, %d edges", vcount(g), ecount(g)), 
       col='darkblue')
  box(col='blue')
  par(mai = mai)
}

findpath <- function(v1, v2){
  
  path <- shortest_paths(mg, from=v1, to=v2)$vpath
  for(i in 1:length(path)){
    p <- path[i]
    for(w in word[cluster.set[[p]]$cluster]) cat(w, ' ')
    cat('\n')
    if(i==length(path)) break
    s1 <- word[cluster.set[[path[i]]]$cluster]
    s2 <- word[cluster.set[[path[i+1]]]$cluster]
    for(w in s1) if(w %in% s2) cat('\n--> ', w, ' -->\n\n')
  } 
}
