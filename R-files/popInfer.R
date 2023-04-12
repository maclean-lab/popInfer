library(glmnet)
source("R-files/Pseudocell-binning.R")
source("R-files/Network-Inference.R")
source("R-files/LASSO-model.R")
source("R-files/Averaged-gene-expression-partition.R")


popInfer <- function(rna.data, gene.access.data, pseudotime, 
                     n.pseudo, alpha, output.directory, printProgress){
  
  if( nrow(rna.data) != nrow(gene.access.data) | ncol(rna.data) != ncol(gene.access.data)) stop('expression matrix and gene accessibility matrix dimensions do not match')
  
  dir.create(output.directory, showWarnings = FALSE)
  
  pseudocell.bins <- PseudocellBinning(pseudotime[,1], pseudotime[,2], 
                                       n.pseudo)
  write.csv(pseudocell.bins, paste0(output.directory, "/pseudocell-bins.csv"))

  avg.expression <- averageCellsPartition(rna.data, pseudocell.bins)
  write.csv(avg.expression, paste0(output.directory, "/pseudocell-expression-matrix.csv"))
  
  avg.accessibility <- averageCellsPartition(gene.access.data, pseudocell.bins)
  write.csv(avg.accessibility, paste0(output.directory, "/pseudocell-gene-accessibility-matrix.csv"))
  
  weight.matrix <- scNetInf(avg.expression, avg.accessibility, alpha, printProgress)
  write.csv(weight.matrix, paste0(output.directory, "/popInfer-weight-matrix.csv"))
}

  
  
  


