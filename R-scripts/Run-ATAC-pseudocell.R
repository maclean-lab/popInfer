library("devtools")
library("BiocManager")
library("ArchR")
library(dineq)
library(data.table)
addArchRGenome("mm10")

source("functions/Janssens-Score-Functions.R")
source("functions/Compute-Janssens-Score.R")

#### load archr object ####
proj <- loadArchRProject(path = "../data/Save-Multiome-Peaks-Strict", force = FALSE)

peakMatrix <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = F,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = NULL
)

## find peakMatrix row index in peakSet 
new_idx_mat <- matrix(data = NA, nrow = nrow(peakMatrix), ncol = 2)
colnames(new_idx_mat) <- c("peak_matrix", "peak_set")
seqnames_vector <- proj@peakSet@seqnames
chr_row_vector <- paste(rep(1, length(table(seqnames_vector)[match(names(table(as.character(seqnames_vector))), names(table(seqnames_vector)))])),
                        table(seqnames_vector)[match(names(table(as.character(seqnames_vector))), names(table(seqnames_vector)))],
                        sep = ":")
chr_row_vector <- unlist(sapply(chr_row_vector, function(x) eval(parse(text=x))))
names(chr_row_vector) <- NULL
chr_row_vector <- paste0(rep(names(table(seqnames_vector)[match(names(table(as.character(seqnames_vector))), names(table(seqnames_vector)))]), table(seqnames_vector)[match(names(table(as.character(seqnames_vector))), names(table(seqnames_vector)))]),
                         ":", chr_row_vector)
new_idx_mat[,1] <- chr_row_vector
new_idx_mat[,2] <- paste0(rep(names(table(seqnames_vector)), table(seqnames_vector)),
                          ":", proj@peakSet$idx)
## GMP branch 
branch.name <- "GMP"

## MEP branch 
#branch.name <- "MEP"

## load pseudocell bin information
pseudocell.bin.df <- read.csv(paste0("outputs/pseudocell-binning-df-", branch.name, ".csv"),
                              row.names = 1)

## define "sample.barcode" to be barcodes that match the ArchR proj format 
pseudocell.bin.df$sample.barcode <- paste0("O_AL#", pseudocell.bin.df$cell.barcodes)

## define and save differentially expressed genes
vec.genes <- read.csv(paste0("../data/", branch.name, "-branch-genes.csv"), row.names = 1)
vec.genes <- vec.genes[,1]

## 5kb upstream from TSS
upstream.extend <- 5000

## output ATAC pseudocells
pseudo.ATAC <- Janssens_GAS(proj, pseudocell.bin.df, vec.genes, upstream.extend, branch.name)



