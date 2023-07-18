#' pcaExtractVariance: Extract variance of each principal component
#' 
#' \code{pcaExtractVariance} constructs, for example, DESeq2 dds object from Salmon count files, condition table and design
#' @param pca_results object return by running prcomp()

# extract variance of each principal component
pcaExtractVariance <- function(pca_results) {
  # adapted from factoextra::get_eig()
  if (inherits(pca_results, "prcomp")) {
    # extract standard deviation of each principal component
    # std_dev <- pca_results$sdev
    
    # compute variances (proportions of variance) - eigenvalues
    # eig_vals <- std_dev^2
    # pr_var <- (pca_results$sdev)^2
    eig_vals <- (pca_results$sdev) ^ 2
    
  } else {
    stop(
      "An object of class : ",
      class(pca_results),
      " can't be handled by the function pca_extract_variance()"
    )
  }
  
  # Proportion of variance explained
  prop_varex <- eig_vals / sum(eig_vals)
  
  # Cumulative proportion of variance explained
  cumsum_varex <- cumsum(prop_varex)
  pca_variances <- data.frame(
    pc = 1:length(eig_vals),
    eig_vals = eig_vals,
    prop_varex = prop_varex,
    cumsum_varex = cumsum_varex
  )
  
  colnames(pca_variances) <-
    c("principal_component",
      "eigenvalues",
      "variance",
      "cummulative_variance")
  
  return(pca_variances)
}

#' pcaPlotVariance: Generate pca variance explained plots
#' 
#' \code{pcaPlotVariance} plot variance shows how many PCs are needed in order to explain e.g. 60% variance (var_expl_needed = 0.6)
#' @param pca_results_extracted results from pcaExtractVariance()
#' @param var_expl_needed % variance explained needed
#' @usage pca_results <- pcaExtractVariance(res_pca)
#' @usage var_expl_needed <- 0.6 # 60% variance explained
#' @usage pca_scree_cumVarExpl <- pcaPlotVariance(pca_results, var_expl_needed = var_expl_needed)

pcaPlotVariance <- function(pca_results_extracted, var_expl_needed = 0.6,
                            npcs = 10, addlabels_scree = TRUE,
                            barfill = "steelblue", barcolor = "steelblue", linecolor = "black") {
  
  # adjusted from factoextra::fviz_eig()
  pca_results_extracted$principal_component <- factor(pca_results_extracted$principal_component)
  # extracting only top n principal components
  pca_results_extracted <- pca_results_extracted[1:min(npcs, nrow(pca_results_extracted)), , drop = FALSE]
  
  # generating scree_plot ----
  scree_plot <- ggplot(pca_results_extracted, aes(x = principal_component, y = variance*100, group = 1)) +
    geom_bar(stat = "identity", fill = barfill, color = barcolor) +
    geom_line(color = linecolor, linetype = "solid") + 
    geom_point(shape = 19, color = linecolor) + 
    labs(title = "Scree plot", x = "Principal components", y = "Percentage of variance explained") + 
    theme_minimal() 
  
  if (addlabels_scree) {
    text_labels <- paste0(round(pca_results_extracted$variance * 100, 1), "%")
    scree_plot <- scree_plot + geom_text(label = text_labels, vjust = -0.4, hjust = 0)
  }
  
  # generating cumsum_plot ----
  # calculate minimum components to explain > 60%; var_expl_needed = 0.6 (60%)
  min_PC_var_needed <- which(pca_results_extracted$cummulative_variance > var_expl_needed)[1]
  min_PC_var_needed_exact <- pca_results_extracted$cummulative_variance[min_PC_var_needed]
  
  cumsum_plot <- ggplot(pca_results_extracted, aes(x = principal_component, y = cummulative_variance*100, group = 1)) +
    geom_line() +
    geom_point(size=3, shape=21, fill=linecolor, color="white", stroke=3) + 
    geom_vline(xintercept = min_PC_var_needed, color = "red", lty=2) +
    geom_label(aes(x=(min_PC_var_needed), y=50, label=paste0(min_PC_var_needed," PCs \n >", var_expl_needed*100,"% variance explained")), 
               color="red") +
    labs(title = "", x = "Principal components", y = "Cumulative proportion of variance explained") + 
    theme_minimal()
  
  # generating final merged plot
  pca_var_plot <- ggpubr::ggarrange(scree_plot, cumsum_plot, ncol = 2, nrow = 1)
  
  return(pca_var_plot) 
}

#' pcaCorrPCs: correlate PCs with covariates
#' 
#' \code{pcaCorrPCs} correlate PCs with covariates, function from pcaExplorer::correlatePCs
#' @param pcaobj     PC_corr <- corrCovarPCs(res_pca, cond_data[cond_interest_varPart], pcs = 1:min_components_needed, pval_exact = TRUE)
#' @param coldata covariates for correlation; e.g. cond_data[cond_interest_varPart]
#' @param pcs how many PCs to correlate
#' @param pval_exact a logical indicating whether an exact p-value should be computed.

pcaCorrPCs <- function(pcaobj, coldata, pcs = 1:4, pval_exact=TRUE) {
  # for ties_error=TRUE use exact=FALSE
  coldataTypes <- sapply(coldata, class)
  x <- pcaobj$x
  
  # check if sample order is the same in pcaobj and coldata
  cat("Identical rownames between pcaobj$x and coldata:", identical(rownames(coldata), rownames(x)),"\n")
  
  res <- matrix(NA, nrow = length(pcs), ncol = ncol(coldata))
  colnames(res) <- colnames(coldata)
  rownames(res) <- paste0("PC_", pcs)
  for (i in 1:ncol(res)) {
    for (j in pcs) {
      
      cat("Correlating PC: ", j, " with covariate: ", colnames(res)[i], "\n")
      
      if (coldataTypes[i] %in% c("factor", "character")) {
        if (length(levels(coldata[, i])) > 1) {
          res[j, i] <- kruskal.test(x[, j], coldata[, i])$p.value
        }
      }
      else {
        if (pval_exact) {
          res[j, i] <- cor.test(x[, j], coldata[, i], method = "spearman")$p.value
        } else {
          res[j, i] <- cor.test(x[, j], coldata[, i], method = "spearman", exact = FALSE)$p.value
        }
      }
    }
  }
  
  return(res)
  
}

#' pcaCorrPCsPlot: generate plot for correlation with PCs
#' 
#' \code{pcaCorrPCsPlot} correlate PCs with covariates, function from pcaExplorer::correlatePCs
#' @param PC_corr object generated with pcaCorrPCs() function
#' [] add option of different plots (bar plot, etc.); see original function for pcaCorrPCs

pcaCorrPCsPlot <- function(PC_corr){
  require("dplyr")
  require(ggpubr)
  #pccorrs <- as.data.frame(PC_corr) %>%
  #  tibble::rownames_to_column(., var = "PC") %>%
  #  tidyr::gather(., covariate, p_val, -PC) %>%
  #  dplyr::mutate(PC = gsub(pattern = "PC_", replacement = "", PC), log10_pval = -log10(p_val)) 
  
  # barplot  
  #pccorrs_plot <- ggbarplot(pccorrs, x = "PC", y = "log10_pval", fill = "covariate", color = "covariate", 
  #                        title = "Significance of the relations between PCs vs covariates",
  #                        xlab = "Principal component", ylab = "-log10(pval)",
  #                       palette = "jco", 
  #                       position = position_dodge()) +
  #geom_hline(yintercept = -log10(0.05), color = "black", linetype="dashed", size=1, alpha=0.8) + 
  #geom_text(aes(max(as.numeric(pccorrs$PC)), -log10(0.05), label = "p-val = 0.05", vjust = -1), color = "black", size = 4)
  
  # lollipop plot
  pccorrs_plot <- as.data.frame(PC_corr) %>%
    tibble::rownames_to_column(., var = "PC") %>%
    tidyr::gather(., covariate, p_val, -PC) %>%
    dplyr::mutate(PC = gsub(pattern = "PC_", replacement = "", PC), 
                  log10_pval = -log10(p_val)) %>% 
    ggpubr::ggdotchart(., x = "PC", y = "log10_pval",
               color = "covariate",                                # Color by groups
               #palette = c("#00AFBB", "#E7B800", "#FC4E07"),Custom color palette
               sorting = "descending",                       # Sort value in descending order
               add = "segment",                             # Add segments from y = 0 to dots
               rotate = FALSE,                                # Rotate vertically
               group = "PC",                                # Order by groups
               dot.size = 8,                                 # Large dot size
               label = round(.$log10_pval, 2),                        # Add mpg values as dot labels
               font.label = list(color = "white", size = 9, 
                                 vjust = 0.5),               # Adjust label parameters
               title = "Significance of the relations between PCs vs covariates",
               xlab = "Principal component", ylab = "-log10(pval)",
               ggtheme = theme_pubr()                        # ggplot2 theme
    ) + geom_hline(yintercept = -log10(0.05), color = "black", linetype="dashed", size=0.5, alpha=0.8) + 
    geom_text(aes(max(as.numeric(PC)), -log10(0.05), label = "p-val = 0.05", vjust = -1), 
              color = "black", size = 4, fontface = "plain", family = "Courier")
  return(pccorrs_plot)
}

#' varPartitionEstimate: ADD description
#' 
varPartitionEstimate <- function(transf_object=NULL, fitform_partVariance=NULL, ntop_genes=1000, ncores=12){
  #TODO:
  # save as a separate file
  # define dependencies
  
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  
  fitform_partVariance <- as.formula(fitform_partVariance)
  
  transf_object_counts <- assay(transf_object)
  geneExpr_partVariance <- rnaSelectTopVarGenes(transf_object_counts, ntop = ntop_genes, type = "var") # top 1000 most variable - for consistency with DESeq2::plotPCA()
  info_partVariance <- as.data.frame(colData(transf_object))
  varPart_fit <- variancePartition::fitExtractVarPartModel(geneExpr_partVariance, fitform_partVariance, info_partVariance)
  varPart_fit_sorted <- variancePartition::sortCols( varPart_fit )
  
  varPart_plot <- variancePartition::plotVarPart(varPart_fit_sorted)
  
  varPart_stats <- varPart_plot$data %>%
    dplyr::select(variable, value) %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(median_varExplained = median(value),
                     mean_varExplained = mean(value),
                     IQR_varExplained = IQR(value),
                     max_varExplained = max(value),
                     min_varExplained = min(value)) %>%
    dplyr::mutate_if(is.numeric, round, 1)
  
  varPart_plot_annot <- varPart_plot + 
    geom_text(data = varPart_stats, aes(x = variable, y = 100, 
                                        label = paste0("mean: ", mean_varExplained)), 
              size = 3, vjust = -0.5)
  
  parallel::stopCluster(cl)   
  
  varPart_results <- list(varPart_fit=varPart_fit, 
                          varPart_stats=varPart_stats, 
                          varPart_plot=varPart_plot, 
                          varPart_plot_annot=varPart_plot_annot)
  return(varPart_results)
}

#' generatePCA: ADD description
#' 
generatePCA <- function(transf_object=NULL, cond_interest_varPart=NULL, color_variable=NULL, shape_variable=NULL, ntop_genes=500){
  #transf_object_counts <- assay(transf_object)
  #ntop_genes=nrow(transf_object_counts) 
  pcaData <- DESeq2::plotPCA(transf_object, intgroup=cond_interest_varPart, returnData=TRUE, ntop=ntop_genes) 
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  ggplot(pcaData, aes(PC1, PC2, color=!!sym(color_variable), shape=!!sym(shape_variable))) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + #+ coord_fixed()
    theme_bw() 
  
}

generatePCA_plus_shape <- function (transf_object = NULL, cond_interest_varPart = NULL, 
                                color_variable = NULL, shape_variable = NULL,   ntop_genes = 500) 
{
  pcaData <- DESeq2::plotPCA(transf_object, intgroup = cond_interest_varPart, 
                             returnData = TRUE, ntop = ntop_genes)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pcaData <- cbind(pcaData, list('Cell_line_label_and_exp_id' = paste(pcaData$cell_line_label, pcaData$experiment)))
  ggplot(pcaData, aes(PC1, PC2, color = !!sym(color_variable), 
                      shape = !!sym(shape_variable))) + geom_point(size = 4) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
    ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
    theme_bw()
}