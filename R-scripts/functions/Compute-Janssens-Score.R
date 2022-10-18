## Input:
## (1) proj = ArchR project 
## (2) vec.genes = genes to compute Janssens score for
## (3) pseudocell.bin.df = output of PseudocellBinning(cell.barcodes, pseudotime.vec, n.pseudo)
## (4) upstream.extend = number of bp to extend upstream of the gene's TSS

## Output:
## Returns and writes matrix (vec.genes x n.pseudo) of Janssens accessibility score


Janssens_GAS <- function(proj, pseudocell.bin.df, vec.genes, upstream.extend, branch.name){
  
  
  pseudocell.bin.df$sample.barcode <- paste0("O_AL#", pseudocell.bin.df$cell.barcodes)
  
  ## subset ArchR project cells
  idxCell <- BiocGenerics::which(proj$cellNames %in% pseudocell.bin.df$sample.barcode)
  cellsSample <- proj$cellNames[idxCell]
  proj <- proj[cellsSample, ]
  
  ## add meta data
  pseudocell.bin.df <- pseudocell.bin.df[match(proj$cellNames, pseudocell.bin.df$sample.barcode),]
  proj$pseudocellBin <- pseudocell.bin.df$pseudocell.bin
  
  ## get bins
  vec.cluster <- unique(pseudocell.bin.df$pseudocell.bin)
  vec.cluster <- vec.cluster[order(vec.cluster)]
  
  mtx_gene_GAS <- matrix(data = NA, nrow = length(vec.genes), ncol = length(vec.cluster)+1)
  
  ## Calculate Janssens Gene Accessibility Score 
  for (i in 1:length(vec.genes)) {
    gene.name <- vec.genes[i]
    gene.idx <- which(proj@geneAnnotation$genes$symbol == gene.name)
    gene.info <- proj@geneAnnotation$genes[gene.idx,]
    if (length(gene.info) == 0) {print("gene not found in annotation");next}
    gene.info_extended <- Extend(gene.info, upstream = upstream.extend, downstream = 0)
    gene.width <- gene.info@ranges@width
    # Find if the extended gene range is overlapping with other genes
    idx_peaks_in_gene.info_extended <- find.peak.regions_HSY(proj, gene.info, upstream.extend)
    overlapping.granges <- proj@peakSet[idx_peaks_in_gene.info_extended,]
    if (length(idx_peaks_in_gene.info_extended) == 0) {print("no peaks in gene region");next}
    # Calculate average counts in each peak # 
    region.access <- accessibility.of.cluster_archr(proj, vec.cluster, idx_peaks_in_gene.info_extended)
    #Calculate distance from each peak to the gene TSS
    dist.to.TSS <- rep(0, length(idx_peaks_in_gene.info_extended))
    for (j in 1:length(overlapping.granges@seqnames)){
      if(as.character(gene.info_extended@strand) == "+") {
        # if start of genomic region upstream of gene, dist = start - tss
        if (overlapping.granges@ranges@start[j] > gene.info@ranges@start){
          dist.to.TSS[j] <- overlapping.granges@ranges@start[j] - gene.info@ranges@start
        } else {
          # else genomic region starts downstream of TSS...
          end.of.region <- overlapping.granges@ranges@start[j] + overlapping.granges@ranges@width[j] - 1
          # if the end of the genomic region doesn't overlap the TSS, dist = TSS - end
          if (end.of.region <  gene.info@ranges@start){
            dist.to.TSS[j] <- gene.info@ranges@start - end.of.region
          } else {
            # else, the genomic region overlaps the TSS
            dist.to.TSS[j] <- 0
          }
        }
      } else {
        # if start of genomic region upstream of gene, dist = start - tss
        if (overlapping.granges@ranges@start[j] > gene.info@ranges@start+gene.info@ranges@width-1){
          dist.to.TSS[j] <- overlapping.granges@ranges@start[j] - gene.info@ranges@start
        } else {
          # else genomic region starts downstream of TSS...
          end.of.region <- overlapping.granges@ranges@start[j] + overlapping.granges@ranges@width[j] - 1
          # if the end of the genomic region doesn't overlap the TSS, dist = TSS - end
          if (end.of.region <  gene.info@ranges@start+gene.info@ranges@width-1){
            dist.to.TSS[j] <- gene.info@ranges@start+gene.info@ranges@width-1 - end.of.region
          } else {
            # else, the genomic region overlaps the TSS
            dist.to.TSS[j] <- 0
          }
        }
      }
    }
    w.d <- exp(-dist.to.TSS/((upstream.extend + gene.width)/2)) + exp(-1)
    gini.scores <- z.Gini(region.access)
    w.b <- exp(gini.scores)
    final.weights <- w.d*w.b
    final.gene.scores <- gene.access.score(region.access, final.weights)
    mtx_gene_GAS[i,] <- c(final.gene.scores, gene.width)
    print(paste0(gene.name, " is done"))
  }
  rownames(mtx_gene_GAS) <- vec.genes
  colnames(mtx_gene_GAS) <- c(vec.cluster, "gene_length")
  
  ## write output file
  write.csv(x = mtx_gene_GAS, file = paste0("outputs/pseudocell-ATAC-", branch.name, ".csv"))
  return(mtx_gene_GAS)
  
}