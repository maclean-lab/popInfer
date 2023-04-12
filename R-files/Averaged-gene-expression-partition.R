## Input:
## expression.matrix = matrix of expression values (genes x cells)
## pseudocell.bin.df = output of PseudocellBinning(cell.barcodes, pseudotime.vec, n.pseudo)
## branch.name = name of the sample/branch

## Output:
## returns and writes matrix (genes x n.pseudo) of pseudocell expression values

averageCellsPartition <- function(expression.matrix, pseudocell.bin.df){
  
  nbins <- max(pseudocell.bin.df$pseudocell.bin)
  pseudo.expr.data <- matrix(0, nrow(expression.matrix), nbins)
  
  for (i in 1: nbins){
    # average over the bin
    indx <- pseudocell.bin.df$pseudocell.bin == i
    num.cells <- nrow(pseudocell.bin.df[indx, ])
    out.vec <- rowSums(expression.matrix[, indx])/num.cells
    pseudo.expr.data[,i] <- out.vec
  }
  
  rownames(pseudo.expr.data) <- rownames(expression.matrix)
  return(pseudo.expr.data)
}

