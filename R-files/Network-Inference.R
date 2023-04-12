## Input:
## (1) sc.data = RNA or ATAC matrix of pseudocells
## (2) alpha.list = sequence of alpha values to average results over
## (3) lag = number of lag values for to use as predictors

## Output:
## weight matrix of signed, directed gene-gene relationships

scNetInf <- function(rna.data, atac.data, alpha.list, printProgress){
  
  equi.dist.pt <- 1:ncol(rna.data)
  degs <- rownames(rna.data)
  
  agg.df <- data.frame(gene=character(),
                       betas.RNA=double(),  
                       target.gene.ATAC=character(),
                       alpha = double()) 
  
  genelist <- degs[1:length(degs)]
  for (geneName in genelist){
    
    if(printProgress == TRUE){
      print(paste0("working on ", geneName, " gene ", 
                   which(genelist == geneName), "/", length(genelist)))
    }
    
    if(var(atac.data[geneName,]) != 0){
      lam.coef.RNA <- OptimLambdaManyAlpha(rna.data, atac.data, geneName, alpha.list)
      
      tmp.df <- lam.coef.RNA
      colnames(tmp.df) <- c("betas.RNA", "gene", "alpha")
      tmp.df$target.gene.ATAC <- rep(geneName, nrow(tmp.df))
      agg.df <- rbind(agg.df, tmp.df)
    }
  }
  
  ## absolute value of coefficient nonzero is inferred TP
  agg.df <- agg.df[abs(agg.df$betas.RNA) > 0, ]
  agg.df <- agg.df[!(agg.df$gene == "(Intercept)"),]
  
  M = matrix(0, length(degs), length(degs))
  rownames(M) <- degs
  colnames(M) <- degs

  for (k in 1:nrow(agg.df)){
    ## add 1 for every time series that identifies
    ## that relationship, subtract 1 if negative relationship
    if (agg.df$betas.RNA[k] > 0){
      M[agg.df$gene[k], agg.df$target.gene.ATAC[k]] <- M[agg.df$gene[k], agg.df$target.gene.ATAC[k]] + 1
    }
    if (agg.df$betas.RNA[k] < 0){
      M[agg.df$gene[k], agg.df$target.gene.ATAC[k]] <- M[agg.df$gene[k], agg.df$target.gene.ATAC[k]] - 1
    }
  }
  M <- M/length(alpha.list)
  return(M)
}


