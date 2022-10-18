#### functions ####
Extend <- function(
  x, upstream = 0, downstream = 0, from.midpoint = FALSE
) {
  if (any(strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- strand(x = x) == "+" | strand(x = x) == "*"
  if (from.midpoint) {
    midpoints <- start(x = x) + (width(x = x) / 2)
    new_start <- midpoints - ifelse(
      test = on_plus, yes = upstream, no = downstream
    )
    new_end <- midpoints + ifelse(
      test = on_plus, yes = downstream, no = upstream
    )
  } else {
    new_start <- start(x = x) - ifelse(
      test = on_plus, yes = upstream, no = downstream
    )
    new_end <- end(x = x) + ifelse(
      test = on_plus, yes = downstream, no = upstream
    )
  }
  ranges(x = x) <- IRanges(start = new_start, end = new_end)
  x <- trim(x = x)
  return(x)
}


find.peak.regions_HSY <- function(archrObject, gene.info, upstream.extend, gene.idx) {
  check.dist <- upstream.extend
  while (check.dist > 0) {
    TSS_checkregion <- Extend(resize(gene.info, 1), upstream = check.dist, downstream = -1)
    
    
    overlapping_genes_idx <- findOverlaps(TSS_checkregion, proj@geneAnnotation$genes[-gene.idx,])
    if (length(overlapping_genes_idx) != 0) {
      overlapping_gene_region <- archrObject@geneAnnotation$genes[-gene.idx,][subjectHits(overlapping_genes_idx),]
      check.dist <- GenomicRanges::distance(gene.info, overlapping_gene_region) - 1
    } else {
      check.dist <- 0
    }
  }
  gene.info_extended <- union(gene.info, TSS_checkregion)
  idx_peaks_in_gene.info_extended <- queryHits(findOverlaps(proj@peakSet, gene.info_extended))
  return(idx_peaks_in_gene.info_extended)
}


accessibility.of.cluster_archr <- function(object, celltype, peak_index){
  # rows are genomic regions, columns are clusters
  
  sum.access.data <- matrix(data = 0, nrow = length(peak_index), ncol = length(celltype))
    
  # for each accessible region we are looking at...
  for (j in 1:length(peak_index)){
    # pulling out correct peakMatrix row
    peakMatrix_row <- which(new_idx_mat[,1] == new_idx_mat[peak_index[j],2])
    # accessibility of each region within each cell
    access.counts <- assays(peakMatrix)$PeakMatrix[peakMatrix_row,]
    for (i in 1:length(celltype)){
      # accessibility of a peak for cells in the cluster
      access.counts.s <- access.counts[names(access.counts) %in% object$cellNames[which(object$pseudocellBin ==celltype[i])]]
      # average over the cells in the bin
      sum.access.data[j,i] <- sum(access.counts.s)/length(access.counts.s)
      
    }
  }
  return(sum.access.data)
  close(pb)
}

## Calculate Gini score and z normalize
z.Gini <- function(region.access){
  gini.scores <- rep(0, nrow(region.access))
  for (i in 1:nrow(region.access)){
    if (!is.na(gini.wtd(region.access[i,]))){gini.scores[i] <- gini.wtd(region.access[i,])}
  }
  # z-standardize 
  if (nrow(region.access) > 1){
    gini.mean <- mean(gini.scores)
    gini.sdev <- sd(gini.scores)
    gini.scores <- (gini.scores - gini.mean)/gini.sdev
  } else {
    # if there is only one peak, set gini score to 0
    gini.scores[1] <- 0
  }
  return(gini.scores)
}

## Compute Janssens gene accessibility score
gene.access.score <- function(region.access, final.weights){ 
  gene.scores <- rep(0, ncol(region.access))
  for (z in 1:ncol(region.access)){
    # normalizing by number of peaks == length(final.weights)
    gene.scores[z] <- sum(region.access[,z] * final.weights)/length(final.weights)
  }
  return(gene.scores)
}
