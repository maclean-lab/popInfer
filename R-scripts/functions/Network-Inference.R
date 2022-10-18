## Input:
## (1) sc.data = RNA or ATAC matrix of pseudocells
## (2) alpha.list = sequence of alpha values to average results over
## (3) lag = number of lag values for to use as predictors

## Output:
## weight matrix of signed, directed gene-gene relationships

scNetInf <- function(sc.data, alpha.list, lag, branch.name){
  
  equi.dist.pt <- 1:ncol(sc.data)
  degs <- rownames(sc.data)
  
  agg.df <- data.frame(gene=character(),
                       betas.RNA=double(),  
                       target.gene=character(),
                       alpha = double()) 
  
  genelist <- degs[1:length(degs)]
  for (geneName in genelist){
    
    lam.coef.RNA <- OptimLambdaManyAlpha(sc.data, geneName, alpha.list, lag)
    
    tmp.df <- lam.coef.RNA
    tmp.df <- tmp.df[c("s1", "lam.gene", "alpha")]
    tmp.df <- aggregate(.~lam.gene+alpha,data=tmp.df,FUN=sum)
    colnames(tmp.df) <- c("gene", "alpha", "betas.RNA")
    tmp.df$target.gene <- rep(geneName, nrow(tmp.df))
    agg.df <- bind_rows(agg.df, tmp.df)
  }
  
  ## absolute value of coefficient nonzero is inferred TP
  agg.df <- agg.df[abs(agg.df$betas.RNA) > 0, ]
  agg.df <- agg.df[!(agg.df$gene == "(Interce"),]
  print(agg.df)
  
  M = matrix(0, length(degs), length(degs))
  rownames(M) <- degs
  colnames(M) <- degs
  print(M)
  for (k in 1:nrow(agg.df)){
    ## add 1 for every time series that identifies
    ## that relationship, subtract 1 if negative relationship
    if (agg.df$betas.RNA[k] > 0){
      M[agg.df$gene[k], agg.df$target.gene[k]] <- M[agg.df$gene[k], agg.df$target.gene[k]] + 1
    }
    if (agg.df$betas.RNA[k] < 0){
      M[agg.df$gene[k], agg.df$target.gene[k]] <- M[agg.df$gene[k], agg.df$target.gene[k]] - 1
    }
  }
  
 write.csv(M, paste0("outputs/weight-matrix-", branch.name, ".csv"))
 return(M)
}


