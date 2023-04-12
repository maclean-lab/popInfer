# popInfer

popInfer is a gene regulatory network inference tool for single-cell multiome data implemented in R. 

### Contents

- `data` contains single-cell gene expression values, single-cell gene accessibility scores, and pseudotime values for five differentiation trajectories across four samples: HSC -> Multipotent transition in yAL, yDR, oAL, oDR, and the HSC -> GMP transition in oAL. 
- `R-files` contains all of the popInfer R functions. 
- `outputs` contains the files returned by popInfer: a gene expression pseudocell matrix, a gene accessibility score pseudocell matrix, and a gene-gene relationship weight matrix. 
- `Run-popInfer.R` is a script to demonstrating how to run popInfer using the samples found in `data`. 

### Paper

### Dependencies
[`glmnet`] for LASSO regression
