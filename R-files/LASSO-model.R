## Input: 
## (1) sc.data = RNA or ATAC matrix of pseudocells
## (2) geneName = target gene name
## (3) alpha.list = sequence of alpha values to average results over
## (4) lag = number of lag values for to use as predictors

## Output: 
## coefficients for regularized linear model for target geneName
## for each alpha in alpha.list


OptimLambdaManyAlpha <- function(rna.data, atac.data, geneName, alpha.list){
  
  agg.df <- data.frame(s1=double(),
                       gene=character(),
                       alpha=double())
  
  ## remove the target gene from predictors
  input.data <- rna.data
  input.data <- input.data[!(rownames(input.data) == geneName),]
  
  ## glmnet fits the model for 200 values of lambda
  fit <- glmnet(t(input.data), t(atac.data[geneName,]), nlambda = 200)
  cvfit <- cv.glmnet(t(input.data), t(atac.data[geneName,]), nlambda = 200)
  
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
      lam.coef$gene <- rownames(coef(fit, s = lambda.trivial.model))
      lam.coef$alpha <- rep(alpha.val, nrow(lam.coef))
      agg.df <- rbind(agg.df, lam.coef)
      
    } else{
      
      num.nonzero.coeffs.standardized <- colSums(betas)/optimal.MSE.nonzero.coeffs
      
      df <- data.frame(lambda.vals, MSE.standardized, num.nonzero.coeffs.standardized)
      
      df$diff <- abs(alpha.val*df$MSE.standardized - (1-alpha.val)*df$num.nonzero.coeffs.standardized)
      ## select lambda that minimizes absolute difference 
      ## alpha controls the tradeoff
      optim.lambda <- df$lambda.vals[which(df$diff == min(df$diff))]
      lambda = optim.lambda
      lam.coef <- data.frame(s1 = matrix(coef(fit, s = lambda)))
      lam.coef$gene <- rownames(coef(fit, s = lambda))
      lam.coef$alpha <- rep(alpha.val, nrow(lam.coef))
      agg.df <- rbind(agg.df, lam.coef)
      
    }

      
  }
  
  return(agg.df)
  
}
