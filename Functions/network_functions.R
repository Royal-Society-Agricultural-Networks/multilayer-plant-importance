#### Functions for matrix and network manipulations #### 

## Convert an edgelist to a matrix
edgelist2Matrix <- function(elist){
  elist$lower <- as.character(elist$lower)
  elist$upper <- as.character(elist$upper)
  rows <- unique(elist$lower)
  cols <- unique(elist$upper)
  network <- matrix(0,nrow = length(rows),ncol = length(cols),dimnames = list(rows,cols))
  for (r in 1:nrow(elist)){
    lowerTaxon <- elist$lower[r]
    upperTaxon <- elist$upper[r]
    network[lowerTaxon,upperTaxon] <- elist$fortotals[r]
  }
  return(network)
}

## Convert a weighted to a binary matrix
binarise <- function(mat){
  m <- mat 
  m[which(m > 0)] <- 1
  return(m)
}