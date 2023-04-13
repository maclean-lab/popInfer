# popInfer

popInfer is a gene regulatory network inference tool for single-cell multiome data implemented in R. 

### Contents

- `data` contains single-cell gene expression values, single-cell gene accessibility scores, and pseudotime values for five differentiation trajectories across four samples: HSC -> Multipotent transition in yAL, yDR, oAL, oDR, and the HSC -> GMP transition in oAL. 
- `R-files` contains all of the popInfer R functions. 
- `outputs` contains the files returned by popInfer: a gene expression pseudocell matrix, a gene accessibility score pseudocell matrix, and a gene-gene relationship weight matrix. 
- `Run-popInfer.R` is a script to demonstrating how to run popInfer using the samples found in `data`. 

### popInfer inputs and parameters
- `rnaData`: gene x cell matrix containing the normalized expression values of the cells/genes to run popInfer on. 
- `geneAccessData`: gene x cell matrix containing the gene accessibility of the cells/genes to run popInfer on (must be the same order of features and cells as `rnaData`). To run popInfer on RNA data only, set `geneAccessData` to be `rnaData`. 
- `pseudotimeData`: a dataframe with two columns containing (1) cell barcodes in the same order as the columns of `rnaData`/`geneAccessData`, and (2) the corresponding pseudotime values for each cell.
- `alpha`: a sequence of alpha values (all in the range [0,1]) which popInfer will be run over. Alpha values closer to zero will produce more sparse networks while alpha values closer to one will produce more dense networks. 
- `outputDir`: path to the directory to write the results of popInfer.
- `printProgress`: boolean variable that designates whether or not to print the gene/index that is being evaluated. 

### popInfer outputs
- `pseudocell-bins.csv`: matrix of three columns, where the first two are the same as the `pseudotimeData` input, and the third column gives the pseudocell bin to which each cell was assigned.
- `pseudocell-expression-matrix.csv`: 
- `pseudocell-gene-accessibility-matrix.csv`: 
- `popInfer-weight-matrix.csv`: gene x gene matrix containing the weights of gene-gene relationships. Rows correspond to regulator genes while columns correspond to target genes. Negative weights indicate an inhibitory relationship. Absolute weight values can be considered to be the "confidence" in a gene pair interaction. 

### Paper

### Dependencies
[`glmnet`] for LASSO regression
