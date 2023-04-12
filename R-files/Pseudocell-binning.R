## Input:
## cell.barcodes = cell barcodes
## pseudotime.vec = vector of pseudotime values per cell in same order as barcodes
## n.pseudo = number of pseudocells to make

## Output:
## dataframe with columns: cell barcodes, pseudotime values, and 
## pseudocell bin 

PseudocellBinning <- function(cell.barcodes, pseudotime.vec, n.pseudo){
  
  n <- length(pseudotime.vec)
  l <- n %% n.pseudo
  
  bin.size <- floor(n/n.pseudo)
  
  numCellsAvg <- rep(bin.size, n.pseudo)
  numCellsAvg[1:l] <- numCellsAvg[1:l] + 1
  
  df <- data.frame(cell.barcodes, pseudotime.vec)
  df$pseudocell.bin <- rep(0, nrow(df))
  
  ## order from smallest to largest pseudotime value
  pt.dists <- order(df$pseudotime.vec)
  
  for (i in 1:n.pseudo){
    # average over the bin
    start.indx <- sum(numCellsAvg[0:(i-1)]) + 1
    end.indx <- sum(numCellsAvg[0:i])
    df[pt.dists[start.indx:end.indx], "pseudocell.bin"] <- i
  }
  
  return(df)
}


