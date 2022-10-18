library(destiny)
library(Seurat)

source("functions/Pseudocell-binning.R")

## load seurat object and set cell type clusters
scd <- readRDS("../data/Rds/scd")
cell.type.clusters <- scd@meta.data$SCT_snn_res.0.25

## take a look at clusters
Idents(scd) <- "SCT_snn_res.0.25"
DimPlot(scd)

## load pseudotime
pseudotime.info <- read.csv("../data/pseudotime-dpt.csv",
                           row.names = 1)

## GMP branch:
branch.name <- "GMP"
branch.clusters <- c(0,1,4)

## MEP branch:
# branch.name <- "MEP"
#branch.clusters <- c(0,1,2,3)


cell.barcodes <- colnames(scd)[cell.type.clusters %in% branch.clusters]
pseudotime.info <- pseudotime.info[rownames(pseudotime.info) %in% cell.barcodes,]

## define number of pseudocells to make
n.pseudo = 100

## results written to outputs directory
pseudocell.bins <- PseudocellBinning(cell.barcodes, pseudotime.info, n.pseudo, branch.name)



