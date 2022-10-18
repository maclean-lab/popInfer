## Input: 
## (1) sc.data = RNA or ATAC matrix of pseudocells
## (2) geneName = target gene name
## (3) alpha.list = sequence of alpha values to average results over
## (4) lag = number of lag values for to use as predictors

## Output: 
## coefficients for regularized linear model for target geneName
## for each alpha in alpha.list


OptimLambdaManyAlpha <- function(sc.data, geneName, alpha.list, lag){
  input.data <- sc.data
  
  agg.df <- data.frame(s1=double(),
                       lam.gene=character(),
                       gene.lag = character(),
                       alpha=double())
  
  input.data <- input.data[!(rownames(input.data) == geneName),]
  
  ## define the row names of the different lags/genes
  lag.rownames <- rep(NA, lag*nrow(input.data))
  for (i in 1:nrow(input.data)){
    start.indx <- lag*(i-1) + 1
    end.indx <- start.indx + lag - 1
    lag.rownames[start.indx:end.indx] <- paste0(rownames(input.data)[i], "_L", c(1:lag))
  }
  
  ## make lagged matrix
  final.input = matrix(0, lag*nrow(input.data), ncol(input.data)-lag)
  
  ## add rownames from above
  rownames(final.input) <- lag.rownames
  for (i in 1:nrow(input.data)){
    for (j in 1:lag){
      genelag <- paste0(rownames(input.data)[i], "_L", j)
      start.indx <- lag - j + 1
      end.indx <- ncol(input.data) - j
      final.input[genelag,] <- t(input.data[rownames(input.data)[i],start.indx:end.indx])
    }
  }
  
  ## remove target gene from the predictors
  final.input <- final.input[!(grepl(geneName, rownames(final.input))),]
  
  ## glmnet fits the model for 200 values of lambda
  fit <- glmnet(t(final.input), t(sc.data[geneName,(lag+1):ncol(input.data)]), nlambda = 200)
  cvfit <- cv.glmnet(t(final.input), t(sc.data[geneName,(lag+1):ncol(input.data)]), nlambda = 200)
  
  lambda.vals <- fit$lambda
  adj.R.2 <- fit$dev.ratio
  
  betas <- as.matrix(fit$beta)
  betas[abs(betas) > 0] <- 1
  
  ## get largetst lambda for which we are estimating with 
  lambda.trivial.model <- max(lambda.vals[colSums(betas) == 0])
  trivial.model.MSE <- cvfit$cvm[cvfit$lambda == lambda.trivial.model]
  MSE.standardized <- cvfit$cvm/(trivial.model.MSE)
  
  ## standardize by number of nonzero coefficients for the lambda
  ## with optimal MSE
  optimal.MSE.nonzero.coeffs <- sum(betas[, cvfit$lambda == cvfit$lambda.min] )
  for (alpha.val in alpha.list){
    ## if the optimal MSE model is the trivial model, no optimization needed
    if (optimal.MSE.nonzero.coeffs == 0){
      #lam.coef <- as.data.frame(coef(fit, s = lambda.trivial.model))
      lam.coef <- data.frame(s1 = matrix(coef(fit, s = lambda.trivial.model)))
      lam.coef$gene.lag <- rownames(coef(fit, s = lambda.trivial.model))
      lam.coef$lam.gene <- str_sub(lam.coef$gene.lag,1,nchar(lam.coef$gene.lag)-3)
      lam.coef$alpha <- rep(alpha.val, nrow(lam.coef))
      agg.df <- bind_rows(agg.df, lam.coef)
      
    } else{
      num.nonzero.coeffs.standardized <- colSums(betas)/optimal.MSE.nonzero.coeffs
      
      df <- data.frame(lambda.vals, MSE.standardized, num.nonzero.coeffs.standardized)
      
      df$diff <- abs(alpha.val*df$MSE.standardized - (1-alpha.val)*df$num.nonzero.coeffs.standardized)
      ## select lambda that minimizes absolute difference 
      ## alpha controls the tradeoff
      optim.lambda <- df$lambda.vals[which(df$diff == min(df$diff))]
      lambda = optim.lambda
      lam.coef <- data.frame(s1 = matrix(coef(fit, s = lambda)))
      lam.coef$gene.lag <- rownames(coef(fit, s = lambda))
      lam.coef$lam.gene <- str_sub(lam.coef$gene.lag,1,nchar(lam.coef$gene.lag)-3)
      lam.coef$alpha <- rep(alpha.val, nrow(lam.coef))
      agg.df <- bind_rows(agg.df, lam.coef)
      
    }
    
  }
  
  return(agg.df)
  
}
