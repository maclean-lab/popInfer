library(glmnet) 
library(dplyr)
source("functions/Network-Inference.R")
source("functions/VAR-model-many-alpha.R")

## GMP branch:
branch.name <- "GMP"
branch.clusters <- c(0,1,4)

## MEP branch:
# branch.name <- "MEP"
#branch.clusters <- c(0,1,2,3)

## load pseudocell matrix
## gene expression
sc.data <- read.csv(paste0("outputs/pseudocell-RNA-", branch.name, ".csv"),
                    row.names = 1)
## gene accessibility
#sc.data <- read.csv(paste0("outputs/pseudocell-ATAC-", branch.name, ".csv"),
#                    row.names = 1)

## sequence of alpha values to run model for
alpha.list <- seq(0.0, 0.1, 0.01)

## number of time lags of predictors to include 
lag = 5

## output weight matrix
weight.matrix <- scNetInf(sc.data, alpha.list, lag, branch.name)

