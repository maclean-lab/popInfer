## Input:
## expression.matrix = matrix of expression values (genes x cells)
## dpt.vec = pseudotime values of cells
## n.pseudo = number of pseudocells to make

## Output:
## matrix (n.pseudo x gene) of pseudocell expression values

averageCellsPartition <- function(expression.matrix, dpt.vec, n.pseudo){
  
  # rows are genes, columns are pseudotime grid points
  pseudo.expr.data <- matrix(0, nrow(expression.matrix), n.pseudo)
  
  n <- length(dpt.vec)
  l <- n %% n.pseudo
    
  bin.size <- floor(n/n.pseudo)
  
  numCellsAvg <- rep(bin.size, n.pseudo)
  numCellsAvg[1:l] <- numCellsAvg[1:l] + 1
  
  ## order from smallest to largest pseudotime value
  pt.dists <- order(dpt.vec)
  
  for (i in 1:n.pseudo){
    # average over the bin
    start.indx <- sum(numCellsAvg[0:(i-1)]) + 1
    end.indx <- sum(numCellsAvg[0:i])
    out.vec <- rowSums(expression.matrix[, pt.dists[start.indx:end.indx]])/numCellsAvg[i]
    pseudo.expr.data[,i] <- out.vec
  }
  
  rownames(pseudo.expr.data) <- rownames(expression.matrix)
  return(pseudo.expr.data)
}

