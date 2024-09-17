#' filterDatasets: filters lowly expressed genes from (DESeq2 object at the moment)
#' 
#' \code{filterDatasets} filters lowly expressed genes either by expression in at least N samples or % of samples
#' @param dds_object DESeq2 object
#' @param abs_filt whether to filter on absolute number of samples.
#' @param abs_filt_samples what is the minimum number of samples for absolute filtering.
#' @param relat_filt fraction of samples for relative filtering.

filterDatasets <- function(dds_object = NULL, 
                           abs_filt = TRUE, 
                           abs_filt_samples = 3, 
                           relat_filt = 0.2) {
    require("DESeq2") # and dependencies to work with dds object
    # filtering samples based on fpm
    # abs_filt = TRUE, filtering on absolute number of samples defined in abs_filt_samples
    # abs_filt = FALSE, filtering on relative number of samples
    
    if (!is.null(dds_object)) {
      cat("Original dds object samples: ", ncol(dds_object), " genes: ", nrow(dds_object), "\n")
      
      if (abs_filt == TRUE) {
        # at least 1 fpm in at least n-number of samples
        # by default considering smallest condition (cluster) with 3 samples
        
        cat("Minimum number of samples with expression:", abs_filt_samples, "\n")
        keep_genes_idx <- rowSums(DESeq2::fpm(dds_object, robust = TRUE) > 1) >= abs_filt_samples
        cat("Number of filtered genes:", sum(keep_genes_idx == FALSE), "\n")
        
        dds_object_filt <- dds_object[keep_genes_idx,]

        cat("Filtered dds object has samples:", ncol(dds_object_filt), "genes:", nrow(dds_object_filt), "\n")
        
        return(dds_object_filt)
        
      } else if (abs_filt == FALSE) {
        # at least 1 fpm in at least x% of samples
        cat(relat_filt * 100, "% of samples:", (relat_filt) * ncol(dds_object), "\n")
        keep_genes_idx <- rowSums(DESeq2::fpm(dds_object, robust = TRUE) > 1) >= relat_filt * ncol(dds_object)
        cat("Number of filtered genes:", sum(keep_genes_idx == FALSE),"\n")

        dds_object_filt <- dds_object[keep_genes_idx,]

        cat("Filtered dds object has samples:", ncol(dds_object_filt), "genes:", nrow(dds_object_filt), "\n")
        
        return(dds_object_filt)
        
      } else {
        cat("abs_filt should be TRUE or FALSE.\n")
      }
      
    } else {
      cat("Need to specify dds_object; e.g. dds_object=dds \n")
    }
  }