#' rnaSelectTopVarGenes: Selects top most variable genes
#' 
#' \code{rnaSelectTopVarGenes} selects top N most variable genes selected by either Median Absolute Deviation (mad) or Variance (var)
#' @param normCounts data.frame with normalized read counts 
#' @param ntop how many most variable genes to select
#'
#' @importFrom genefilter rowVars rowSds

# extracting using mad
#TODO:
# [ ] -  add option to use rowVars as in DESeq2:::plotPCA.DESeqTransform
rnaSelectTopVarGenes <- function(normCounts=NULL, 
                                 ntop=500, 
                                 type="mad", 
                                 mad_cutoff=NULL){
  #TODO:
  # [ ] - add option to plot mad as a function of ordered number of genes
  # [ ] add option to select based on variance
  # [ ] add functions without dependencies; re-implement from genefilter!
  # [ ] unify mad and var approaches
  # [ ] add rowSds! 
  
  if (is.null(normCounts)) {
    stop("Need to provide normCounts matrix!")
  } 
  else {
    if (type=="mad") {
      
      row_mads <- apply(normCounts, 1, mad) # top variable genes, measured by mad (median absolute deviation).
      
      if (!is.null(mad_cutoff)) {
        cat("Selecting Ntop genes based on mad_cutoff: ", mad_cutoff, "\n")
        ords <- names(sort(row_mads[row_mads > mad_cutoff], decreasing=TRUE)) #[row_mads > mad_cutoff]
        select <- normCounts[ords, ]
      } else {
        cat("Selecting", ntop, "Ntop genes based on mad.\n")
        ords <- names(sort(row_mads, decreasing=TRUE))
        select <- normCounts[ords[1:ntop], ]
      }
      
    } else if (type=="var"){
      cat("Selecting", ntop, "Ntop genes based on variance.\n")
      
      #row_vars <- apply(normCounts, 1, var) # top variable genes, measured by variance.
      row_vars <- genefilter::rowVars(normCounts)
      
      order_select_vars <- order(row_vars, decreasing = TRUE)[seq_len(min(ntop, length(row_vars)))]
      select <- normCounts[order_select_vars, ]
      
    } else {
      stop("Need to select type='mad' or type='var'\n")
    }
    
    return(select)
    
  }
}