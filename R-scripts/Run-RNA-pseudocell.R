library(destiny)
library(Seurat)

source("functions/Averaged-gene-expression-partition.R")

## load seurat object and set cell type clusters
scd <- readRDS("../data/Rds/scd")
cell.type.clusters <- scd@meta.data$SCT_snn_res.0.25

## GMP branch:
branch.name <- "GMP"
branch.clusters <- c(0,1,4)

## MEP branch:
# branch.name <- "MEP"
#branch.clusters <- c(0,1,2,3)

## load pseudocell bin information
pseudocell.bin.df <- read.csv(paste0("outputs/pseudocell-binning-df-", branch.name, ".csv"),
                              row.names = 1)

## define and save differentially expressed genes
degs <- read.csv(paste0("../data/", branch.name, "-branch-genes.csv"), row.names = 1)
degs <- degs[,1]

## define branch expression matrix
expression.matrix <- as.matrix(scd@assays$SCT@data)[, cell.type.clusters %in% branch.clusters]
expression.matrix <- expression.matrix[rownames(expression.matrix) %in% degs,]

## output pseudocell RNA expression 
avg.expression <- averageCellsPartition(expression.matrix, pseudocell.bin.df, branch.name)

