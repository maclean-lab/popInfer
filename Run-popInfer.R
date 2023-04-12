## dependencies: glmnet
source("R-files/popInfer.R")

## specify dataset 
## HSC -> Multipotent transition (yAl, yDR, oAL, oDR)
## HSC -> GMP transition (oAL-GMP)
sampleName <- "yAL"

## specify pseudotime cell cycle status 
## ("cell-cycle" or "cell-cycle-regressed")
cellCycleStatus <- "cell-cycle"

## specify output directory name
outputDir <- paste0("outputs/output-", sampleName, "-", cellCycleStatus)

## specify number of pseudocells
numberPseudocells <- 80

## specify alpha sequence
alpha <- seq(0.0, 0.4, 0.001)

## load RNA matrix
rnaData <- read.csv(paste0("data/", sampleName, "/gene-expr-data-", 
                            sampleName,".csv"), row.names = 1)

## load gene accessibility score matrix
geneAccessData <- read.csv(paste0("data/", sampleName, "/gene-access-scores-", 
                                    sampleName,".csv"), row.names = 1)

## load pseudotime and barcode data
pseudotimeData <- read.csv(paste0("data/", sampleName, "/",
                                  cellCycleStatus, "-pseudotime-", sampleName, 
                                  ".csv"), row.names = 1)

## run popInfer
popInfer(rnaData, geneAccessData, pseudotimeData, 
         numberPseudocells, alpha, outputDir, printProgress = TRUE)




