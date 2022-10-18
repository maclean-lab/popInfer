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

## GMP branch (clusters = 0,1,4):
branch.name <- "GMP"

## MEP branch (clusters = 0,1,2,3):
#branch.name <- "MEP"

cell.barcodes <- read.csv(paste0("../data/", branch.name, "-branch-cells.csv"), row.names = 1)
cell.barcodes <- cell.barcodes$x

pseudotime.info <- pseudotime.info[rownames(pseudotime.info) %in% cell.barcodes,]

## define number of pseudocells to make
n.pseudo = 100

## results written to outputs directory
pseudocell.bins <- PseudocellBinning(cell.barcodes, pseudotime.info, n.pseudo, branch.name)



