# Results function
# components:
# 1. extract comparison, numerator and denominator
# 2. supp. func: collect mean norm counts for each numerator and denominator
# 3. 
# [ ] use apeglm, or ashr (for contrast) instead of normal shrinkage and plot comparisons
# [ ] specify coefficient  number rather than name
# https://www.bioconductor.org/packages/devel/bioc/vignettes/apeglm/inst/doc/apeglm.html
# Be aware that DESeq2’s lfcShrink interface provides LFCs on the log2 scale, while apeglm provides coefficients on the natural log scale.
# [ ] add MA-plot!
# [ ] add genomic ranges for visualiztion with Gviz
# generateResults: add imported functions!


#' extractResults
#'
#' @param dds_filt_extractResults 
#' @param ensemblAnnot 
#' @param padj_cutoff 
#' @param log2FC_cutoff 
#' @param results_name 
#' @param OUTPUT_DIR 
#' @param condition_test 
#' @param variable 
#' @param shrinkage_type
#'
#' @return
#' @export
#'
#' @examples
#' results_A549 <- extractResults(dds_filt_extractResults = dds_filt_A549, 
#'                                ensemblAnnot = ensemblAnnot, 
#'                                padj_cutoff = padj_cutoff, 
#'                                log2FC_cutoff = log2(1.3), 
#'                                results_name = "A549", 
#'                                OUTPUT_DIR=OUTPUT_DIR, 
#'                                condition_test="condition_siRIC_vs_ctl_ntRIC", 
#'                                variable="condition")
#'                                
extractResults <- function(dds_filt_extractResults=NULL,
                           ensemblAnnot=NULL,
                           padj_cutoff=NULL,
                           log2FC_cutoff=NULL,
                           results_name = NULL,
                           OUTPUT_DIR=NULL,
                           condition_test=NULL, 
                           variable=NULL,
                           shrinkage_type="normal",
                           use_contrast=FALSE){
  # add require(DESeq2)
  
  
  # condition_test <- "condition_siRIC_vs_ctl_ntRIC"
  # variable <- "condition"
  # results_name <- cell_line # something that distinguishes results from other; used in dir_name
  
  # creating unique results DRI
  RESULTS_DIR = paste0(OUTPUT_DIR, results_name, "/")
  cat("Saving results in RESULTS_DIR:", RESULTS_DIR, "\n")
  dir.create(RESULTS_DIR, showWarnings = TRUE)
  
  # extracting numerator/denominator 
  cond_numerator <- gsub(pattern = paste0("(", variable,"_)(.+)(_vs_.+)"), replacement = "\\2", condition_test)
  cond_denominator <- gsub(pattern = paste0("(", variable,"_)(.+_vs_)(.+)"), replacement = "\\3", condition_test)
  
  # extracting results and annotating
  results_all <- generateResults(dds_object=dds_filt_extractResults, 
                                 res_extract=condition_test, 
                                 annot_df=ensemblAnnot,
                                 variable=variable, 
                                 padj_cutoff=padj_cutoff, 
                                 shrinkage_type=shrinkage_type,
                                 use_contrast=use_contrast) # "apeglm", "normal"
  results_data <- results_all$results_data
  results_annot <- results_all$results_data_annot
  
  # extracting signigicant results
  results_annot_signif <- results_annot %>%
    dplyr::filter((!is.na(padj) & (padj < padj_cutoff)) & abs(log2FoldChange) > log2FC_cutoff)
  
  # for metacore
  #results_annot_signif_geneGO <- results_annot_signif %>%
  #  dplyr::select(ensembl_id, log2FoldChange, padj)
  
  # saving results for use in notebook
  cat("Saving results...\n")
  cat(paste0(cond_numerator,"_vs_", cond_denominator,"_dds_filt_extractResults.RData\n"))
  save(dds_filt_extractResults, results_data, results_annot,
       file = paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_dds_filt_extractResults.RData"))
  #writeLines(capture.output(summary(results_data, alpha = padj_cutoff)), paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_summary.csv"))
  cat(paste0(cond_numerator,"_vs_", cond_denominator,".(csv)(xlsx)\n"))
  readr::write_csv(results_annot, path=paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,".csv"))
  writexl::write_xlsx(results_annot, path=paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,".xlsx"))
  cat(paste0(cond_numerator,"_vs_", cond_denominator,"_signif.(csv)(xlsx)\n"))
  readr::write_csv(results_annot_signif, path=paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_signif.csv"))
  #write_tsv(results_annot_signif_geneGO, path=paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_signifGeneGO.tsv"))
  writexl::write_xlsx(results_annot_signif, path=paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_signif.xlsx"))
  #write_tsv(results_annot_signif_geneGO, path=paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_signifGeneGO.tsv"))
  
  # ?! using result_names in files
  # save(dds_filt_extractResults, results_data, results_annot, 
  #      file = paste0(RESULTS_DIR, results_name, "_", cond_numerator,"_vs_", cond_denominator,"_dds_filt_extractResults.RData"))
  # #writeLines(capture.output(summary(results_data, alpha = padj_cutoff)), paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_summary.csv"))
  # cat(paste0(cond_numerator,"_vs_", cond_denominator,".(csv)(xlsx)\n"))
  # readr::write_csv(results_annot, path=paste0(RESULTS_DIR, results_name, "_", cond_numerator,"_vs_", cond_denominator,".csv"))
  # writexl::write_xlsx(results_annot, path=paste0(RESULTS_DIR, results_name, "_", cond_numerator,"_vs_", cond_denominator,".xlsx"))
  # cat(paste0(cond_numerator,"_vs_", cond_denominator,"_signif.(csv)(xlsx)\n"))
  # readr::write_csv(results_annot_signif, path=paste0(RESULTS_DIR, results_name, "_", cond_numerator,"_vs_", cond_denominator,"_signif.csv"))
  # #write_tsv(results_annot_signif_geneGO, path=paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_signifGeneGO.tsv"))
  # writexl::write_xlsx(results_annot_signif, path=paste0(RESULTS_DIR, results_name, "_", cond_numerator,"_vs_", cond_denominator,"_signif.xlsx"))
  # #write_tsv(results_annot_signif_geneGO, path=paste0(RESULTS_DIR, cond_numerator,"_vs_", cond_denominator,"_signifGeneGO.tsv"))
  
  return(list(results_annot = results_annot,
              results_annot_signif = results_annot_signif))  
  
}


#' generateResults
#'
#' @description 
#' Generates results data.frame and annotated results data.frame
#'
#' @param dds_object DESeq2::DESeqDataSet object
#' @param res_extract DESeq2::resultsNames 
#' @param annot_df annotation data.frame as generated by generateEnsemblAnnotation
#' @param OUTPUT_DIR output directory - NOT IMPLEMENTED AT THE MOMENT
#' @param variable variable of interest - e.g. condition, treatment; complex designs NOT IMPLEMENTED AT THE MOMENT
#' @param padj_cutoff alpha in DESeq2::results function; the significance cutoff used for optimizing the independent filtering
#' @param shrinkage_type shrinkage type used in DESeq2::lfcShrink
#'
#' @return
#' results_list that contains results data.frame and annotated results data.frame
#' @export
#'
#' @examples
#' results_all <- generateResults(dds_object=dds_filt_extractResults, 
#'                                res_extract="condition_siRIC_vs_ctl_ntRIC", 
#'                                annot_df=ensemblAnnot,
#'                                variable="condition", 
#'                                padj_cutoff=0.05, 
#'                                shrinkage_type="normal") # "apeglm", "normal"
#'                            
generateResults <- function(dds_object=NULL, 
                            res_extract=NULL, 
                            annot_df=ensemblAnnot, 
                            OUTPUT_DIR=NULL, 
                            variable="condition", 
                            padj_cutoff=0.05, 
                            shrinkage_type="apeglm",
                            use_contrast=FALSE){
  require("DESeq2")
  require("dplyr")
  require("tibble")
  
  res_extract <- res_extract
  # res_extract - what results to extract; e.g. res_extract = results_coeff[[res_test]]
  # OUTPUT_DIR - output dir for the experiment; not needed at the moment!
  # [ ] add validity and object checks
  # shrinkage_type - "normal", "apeglm", "ashr"; currently "apeglm" doesn't work with contrast need to revert back to "normal" (or use "ashr")
  # annot_df - data.frame that contains annotation as generated during preProcessing step; annot_df=ensemblAnnot
  
  # padding stringr::str_pad(stringr::str_pad("Extracting Results", width = 20, side="both"), 50, side="both", pad="-")
  cat(paste(rep("-", 20), collapse = ""), "Extracting Results", paste(rep("-", (30-nchar("Extracting Results"))), collapse = ""), "\n")
  cat("Analysing:", res_extract,"\n") 
  
  cond_numerator <- gsub(pattern = paste0("(", variable,"_)(.+)(_vs_.+)"), replacement = "\\2", res_extract)
  cond_denominator <- gsub(pattern = paste0("(", variable,"_)(.+_vs_)(.+)"), replacement = "\\3", res_extract)
  
  cat("Comparison:",cond_numerator, "vs", cond_denominator, "\n")
  
  # extracting normalized counts for conditions tested and adding mean expression per group
  normalized_counts_AddedMean <- meanExprsPerGroup(dds_object=dds_object, 
                                                   cond_numerator=cond_numerator, 
                                                   cond_denominator=cond_denominator,
                                                   variable = variable)
  
  # Extracting results ----
  # name and contrast should be equivalent, except for (http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#theory): 
  # contrast will additionally set to 0 the estimated LFC in a comparison of two groups, where all of the counts in the two groups are equal to 0 (while other groups have positive counts).
  
  
  if (use_contrast | !(res_extract %in% resultsNames(dds_object))) {
    cat("using contrast to extract results\n")
    # use contrast when results not in resultsNames or if use of contrast is specified by use_contrast=TRUE
    # if apeglm - redefine base factor so apeglm can be used with coef
    res_table_unshrunken <- DESeq2::results(dds_object, 
                                            contrast=c(variable, cond_numerator, cond_denominator), 
                                            parallel = TRUE, alpha = padj_cutoff)
    if(shrinkage_type == "apeglm"){
      cat("apeglm is not implemented for contrast, using: normal or ashr (or specify shrinkage_type='normal') instead\n")
      shrinkage_type <- "normal" # or specify ashr instead
      res_table <- DESeq2::lfcShrink(dds_object, 
                                     contrast=c(variable, cond_numerator, cond_denominator), 
                                     res=res_table_unshrunken, type = shrinkage_type)
    } else {
      res_table <- DESeq2::lfcShrink(dds_object, 
                                     contrast=c(variable, cond_numerator, cond_denominator), 
                                     res=res_table_unshrunken, type = shrinkage_type)
    }
  } else if (res_extract %in% resultsNames(dds_object)) {
    cat("name exists in resultsNames(dds_object), using name to extract results\n")
    res_table_unshrunken <- DESeq2::results(dds_object, name=res_extract, alpha = padj_cutoff)
    cat("Adding shrunken log2 fold changes (LFC) and SE to a results table...using:", shrinkage_type,"\n")
    res_table <- DESeq2::lfcShrink(dds_object, coef=res_extract, res=res_table_unshrunken, type = shrinkage_type)
    
  } else {
    # add more descriptive error message
    stop("There is a problem with results extraction!")
  }
  
  # Constructing final results to be exported ----
  # calculate FC from log2FC, also add other columns
  # !!! when running apeglm column "stat" is not generated 
  
  results_data_annot <- as.data.frame(res_table) %>%
    tibble::rownames_to_column("ensembl_id") %>% 
    dplyr::mutate(FoldChange=ifelse(log2FoldChange < 0, -2^(abs(log2FoldChange)), 2^(log2FoldChange))) %>% #adding FoldChange column
    dplyr::left_join(annot_df, "ensembl_id") %>% # adding annotationl entrez_ids, gene symbols,...
    dplyr::left_join(normalized_counts_AddedMean, "ensembl_id") %>% # adding normalized counts and average counts per group
    dplyr::arrange(padj) %>% # ordering based on p.adj
    dplyr::select(ensembl_id, contains("_symbol"), contains("description"), entrezgene,
                  baseMean, MeanExpr_denominator, MeanExpr_numerator, log2FoldChange, lfcSE, FoldChange, contains("stat"), pvalue, padj, 
                  gene_biotype, chromosome_name,  start_position, end_position, strand, 
                  starts_with(paste0(cond_denominator,"-")), starts_with(paste0(cond_numerator,"-"))) %>% # everything() - other columns not mentioned; starts_with(paste0(cond_denominator,"_S") may not work if other references names dont follow with _S
    dplyr::rename_at(., .vars = "MeanExpr_numerator", .funs = funs(gsub("numerator", "", paste0("MeanExpr_", cond_numerator)))) %>% # find better way of renaming!  
    dplyr::rename_at(., .vars = "MeanExpr_denominator", .funs = funs(gsub("denominator", "", paste0("MeanExpr_", cond_denominator))))
  
  cat(paste(rep("-", 20), collapse = ""), "Finished", paste(rep("-", (30-nchar("Finished"))), collapse = ""), "\n")
  #cat(paste(rep("-", 20), collapse = ""), "Finished", paste(rep("-", 20), collapse = ""), "\n")    
  
  results_list = list(results_data = res_table, 
                      results_data_annot = results_data_annot)
  
  return(results_list)
  
}

meanExprsPerGroup <- function(dds_object=NULL, 
                              #condition_test=NULL,
                              variable=NULL,
                              cond_numerator=NULL, 
                              cond_denominator=NULL){
  require(DESeq2)
  require(dplyr)
  # extracts expression matrix (all samples) for numerator and denominator
  # calculates mean expression per group (numerator, denominator)
  # [ ] add check that there is filenames and condition column!
  # Normalized counts and means per group
  # creating only subset for actual comparison
  # may need to change condition to some other variable
  # res_extract <- condition_test
  # cond_numerator <- gsub(pattern = paste0("(", variable,"_)(.+)(_vs_.+)"), replacement = "\\2", res_extract)
  # cond_denominator <- gsub(pattern = paste0("(", variable,"_)(.+_vs_)(.+)"), replacement = "\\3", res_extract)
  
  new_sample_names <- as.data.frame(colData(dds_object)) 
  
  if("filenames" %in% names(new_sample_names)){
    new_sample_names <- new_sample_names %>%
      dplyr::select(filenames, tidyselect::all_of(variable)) 
    
  } else {
    new_sample_names <- new_sample_names %>%
      tibble::rownames_to_column(., var = "filenames") %>%
      dplyr::select(filenames, tidyselect::all_of(variable)) 
  }
  
  # check if variable is factor otherwise arranging in the next step will be "random"
  #if (is.factor()) {
  #  
  #}
  
  # extract filenames for each of the conditions and pivot_table
  new_sample_names <- new_sample_names %>%
    dplyr::select(filenames, !!as.name(variable)) %>% # selecting filenames and variable of interest (e.g.)
    dplyr::filter(!!as.name(variable) == cond_numerator | !!as.name(variable) == cond_denominator) %>% # filtering to keep only numerator and denominator samples
    dplyr::transmute(filenames,
                     denom_num_extract = factor(!!as.name(variable), levels=c(cond_denominator, cond_numerator))) %>% # re-factoring denominator, then numerator (but this should have been done in cond_data already)
    dplyr::arrange(., denom_num_extract) %>% # convert the strings to names with as.name and !! unquote (bang-bang) !!as.name(variable)
    # reorder variable according to numerator, denominator; so the expression output is in correct order - this assumes previous correct ordering
    dplyr::mutate(new_name = paste(denom_num_extract, filenames, sep="-")) # creating new column with variable of interest and filename
  
  
  normalized_counts <- NULL # just to make sure it does not exist from previous run
  normalized_counts <- data.frame(counts(dds_object, normalized = TRUE))
  # extracting subset
  normalized_counts <- normalized_counts %>%
    dplyr::select(new_sample_names$filenames) # keeping only conditions that are being compared and following previous order
  
  colnames(normalized_counts) <- new_sample_names$new_name
  normalized_counts <- normalized_counts %>%
    tibble::rownames_to_column("ensembl_id")
  
  normalized_counts_AddedMean <- normalized_counts %>%
    dplyr::mutate(., 
                  MeanExpr_denominator = rowMeans(dplyr::select(., matches(paste0(cond_denominator,"-"))), na.rm = TRUE),
                  MeanExpr_numerator = rowMeans(dplyr::select(., matches(paste0(cond_numerator,"-"))), na.rm = TRUE)) # more robust regex?!
  
  # rename mean_numerator, mean_denominator in the final column
  # use rename_all! rename(new_sample_names$new_sample_names)
  # rename colnames to condition_Sample number name
  # calculate mean across normalized counts
  
  return(normalized_counts_AddedMean)
}




test_lfcThreshold=FALSE
if(test_lfcThreshold){
  
  AAA <- DESeq2::results(dds_filt_extractResults, name=results_coeff[[res_test]], lfcThreshold=0.58, alpha = padj_cutoff)
  AAA_df <- as.data.frame(AAA)
  BBB <- DESeq2::results(dds_filt_extractResults, name=results_coeff[[res_test]], lfcThreshold=0, alpha = padj_cutoff)
  BBB_df <- as.data.frame(BBB)
  sum(!is.na(AAA_df$padj) & AAA_df$padj < 0.05)
  sum((!is.na(BBB_df$padj) & BBB_df$padj < 0.05) & (abs(BBB_df$log2FoldChange) > 0.58))
  
  #design(dds_filt) <- ~ condition + library_preparation
  #dds_filt <- DESeq(dds_filt)
  
  AAA_df_filt <- AAA_df %>%
    filter(!is.na(padj) & padj != 1)
  AAA_df_pval_hist <- ggplot(AAA_df_filt, aes(x = pvalue)) +
    #theme_minimal() +
    geom_histogram(binwidth = 0.025, boundary = 0 ) +
    ggtitle("lfcThreshold=0.58")
  
  BBB_df_filt <- BBB_df %>%
    filter(!is.na(padj) & padj != 1)
  BBB_df_pval_hist <- ggplot(BBB_df_filt, aes(x = pvalue)) +
    theme_minimal() +
    geom_histogram(binwidth = 0.025, boundary = 0 ) +
    ggtitle("lfcThreshold=0")
  
  plot(AAA_df$log2FoldChange, BBB_df$log2FoldChange)
  abline(0,1)
}

generateResults_upd <- function(dds_object=NULL, coeff_name=NULL, cond_numerator=NULL,cond_denominator=NULL, padj_cutoff = 0.05, log2FC_cutoff = 0.58, cond_variable=NULL){
  res_table_unshrunken <- DESeq2::results(dds_object, 
                                          name=coeff_name,
                                          parallel = TRUE, alpha = padj_cutoff)
  #renv::install("bioc::apeglm")
  res_table <- DESeq2::lfcShrink(dds_object, 
                                 coef=coeff_name,
                                 res=res_table_unshrunken, type = "apeglm")
  #FIXME:
  # [ ] - fix mean expression part
  normalized_counts_AddedMean <- meanExprsPerGroup(dds_object=dds_object,
                                                   cond_numerator=cond_numerator,
                                                   cond_denominator=cond_denominator,
                                                   variable = cond_variable)
  
  # dplyr::select(ensembl_id, contains("_symbol"), contains("description"), entrezgene,
  #               baseMean, MeanExpr_denominator, MeanExpr_numerator, log2FoldChange, lfcSE, FoldChange, contains("stat"), pvalue, padj, 
  #               gene_biotype, chromosome_name,  start_position, end_position, strand, 
  #               starts_with(paste0(cond_denominator,"-")), starts_with(paste0(cond_numerator,"-"))) #%>% # everything() - other columns not mentioned; starts_with(paste0(cond_denominator,"_S") may not work if other references names dont follow with _S
  # 
  results_data_annot <- as.data.frame(res_table) %>%
    tibble::rownames_to_column("ensembl_id") %>% 
    dplyr::mutate(FoldChange=ifelse(log2FoldChange < 0, -2^(abs(log2FoldChange)), 2^(log2FoldChange))) %>% #adding FoldChange column
    dplyr::left_join(ensemblAnnot, "ensembl_id") %>% # adding annotationl entrez_ids, gene symbols,...
    dplyr::left_join(normalized_counts_AddedMean, "ensembl_id") %>% # adding normalized counts and average counts per group
    dplyr::arrange(padj) %>% # ordering based on p.adj
    dplyr::select(ensembl_id, contains("_symbol"), contains("description"), entrezgene,
                  baseMean, MeanExpr_denominator, MeanExpr_numerator, log2FoldChange, lfcSE, FoldChange, contains("stat"), pvalue, padj,
                  gene_biotype, chromosome_name,  start_position, end_position, strand,
                  starts_with(paste0(cond_denominator,"-")), starts_with(paste0(cond_numerator,"-"))) %>% # everything() - other columns not mentioned; starts_with(paste0(cond_denominator,"_S") may not work if other references names dont follow with _S
    #%>% # everything() - other columns not mentioned; starts_with(paste0(cond_denominator,"_S") may not work if other references names dont follow with _S
    dplyr::rename_at(., .vars = "MeanExpr_numerator", .funs = funs(gsub("numerator", "", paste0("MeanExpr_", cond_numerator)))) %>% # find better way of renaming!  
    dplyr::rename_at(., .vars = "MeanExpr_denominator", .funs = funs(gsub("denominator", "", paste0("MeanExpr_", cond_denominator))))
  
  results_data_annot_signif <- results_data_annot %>%
    dplyr::filter((!is.na(padj) & (padj < padj_cutoff)) & abs(log2FoldChange) > log2FC_cutoff)
  
  temp_results_summary_df <- data.frame(test = coeff_name,
                                        design = paste(as.character(design(dds_object)), collapse=""),
                                        signif_genes = nrow(results_data_annot_signif),
                                        signif_genes_UP = sum(results_data_annot_signif$log2FoldChange > log2FC_cutoff),
                                        signif_genes_DOWN = sum(results_data_annot_signif$log2FoldChange < log2FC_cutoff),
                                        cutoffs = paste0("(!is.na(padj) & (padj < ", padj_cutoff,")) & abs(log2FoldChange) > ", log2FC_cutoff))
  
  
  results_list<-list(results_signif=results_data_annot_signif, de_details=temp_results_summary_df, results_all=results_data_annot)
  return(results_list)
}






generateResults_upd2 <- function(dds_object=NULL, coeff_name=NULL, cond_numerator=NULL,cond_denominator=NULL, padj_cutoff = 0.05, log2FC_cutoff = 0.58, cond_variable=NULL){
  res_table_unshrunken <- DESeq2::results(dds_object, 
                                          contrast=c(cond_variable, cond_numerator, cond_denominator),
                                          parallel = TRUE, alpha = padj_cutoff)
  #renv::install("bioc::apeglm")
  res_table <- DESeq2::lfcShrink(dds_object, 
                                 coef=coeff_name,
                                 res=res_table_unshrunken, type = "apeglm")
  #FIXME:
  # [ ] - fix mean expression part
  normalized_counts_AddedMean <- meanExprsPerGroup(dds_object=dds_object,
                                                   cond_numerator=cond_numerator,
                                                   cond_denominator=cond_denominator,
                                                   variable = cond_variable)
  
  # dplyr::select(ensembl_id, contains("_symbol"), contains("description"), entrezgene,
  #               baseMean, MeanExpr_denominator, MeanExpr_numerator, log2FoldChange, lfcSE, FoldChange, contains("stat"), pvalue, padj, 
  #               gene_biotype, chromosome_name,  start_position, end_position, strand, 
  #               starts_with(paste0(cond_denominator,"-")), starts_with(paste0(cond_numerator,"-"))) #%>% # everything() - other columns not mentioned; starts_with(paste0(cond_denominator,"_S") may not work if other references names dont follow with _S
  # 
  results_data_annot <- as.data.frame(res_table) %>%
    tibble::rownames_to_column("ensembl_id") %>% 
    dplyr::mutate(FoldChange=ifelse(log2FoldChange < 0, -2^(abs(log2FoldChange)), 2^(log2FoldChange))) %>% #adding FoldChange column
    dplyr::left_join(ensemblAnnot, "ensembl_id") %>% # adding annotationl entrez_ids, gene symbols,...
    dplyr::left_join(normalized_counts_AddedMean, "ensembl_id") %>% # adding normalized counts and average counts per group
    dplyr::arrange(padj) %>% # ordering based on p.adj
    dplyr::select(ensembl_id, contains("_symbol"), contains("description"), entrezgene,
                  baseMean, MeanExpr_denominator, MeanExpr_numerator, log2FoldChange, lfcSE, FoldChange, contains("stat"), pvalue, padj,
                  gene_biotype, chromosome_name,  start_position, end_position, strand,
                  starts_with(paste0(cond_denominator,"-")), starts_with(paste0(cond_numerator,"-"))) %>% # everything() - other columns not mentioned; starts_with(paste0(cond_denominator,"_S") may not work if other references names dont follow with _S
    #%>% # everything() - other columns not mentioned; starts_with(paste0(cond_denominator,"_S") may not work if other references names dont follow with _S
    dplyr::rename_at(., .vars = "MeanExpr_numerator", .funs = funs(gsub("numerator", "", paste0("MeanExpr_", cond_numerator)))) %>% # find better way of renaming!  
    dplyr::rename_at(., .vars = "MeanExpr_denominator", .funs = funs(gsub("denominator", "", paste0("MeanExpr_", cond_denominator))))
  
  results_data_annot_signif <- results_data_annot %>%
    dplyr::filter((!is.na(padj) & (padj < padj_cutoff)) & abs(log2FoldChange) > log2FC_cutoff)
  
  temp_results_summary_df <- data.frame(test = coeff_name,
                                        design = paste(as.character(design(dds_object)), collapse=""),
                                        signif_genes = nrow(results_data_annot_signif),
                                        signif_genes_UP = sum(results_data_annot_signif$log2FoldChange > log2FC_cutoff),
                                        signif_genes_DOWN = sum(results_data_annot_signif$log2FoldChange < log2FC_cutoff),
                                        cutoffs = paste0("(!is.na(padj) & (padj < ", padj_cutoff,")) & abs(log2FoldChange) > ", log2FC_cutoff))
  
  
  results_list<-list(results_signif=results_data_annot_signif, de_details=temp_results_summary_df, results_all=results_data_annot)
  return(results_list)
}


generateResults_upd_for_review <- function(dds_object = NULL,
                                           contrast_list = NULL, 
                                           coeff_name = NULL, 
                                           padj_cutoff = 0.05, 
                                           log2FC_cutoff = 0.58){
  #dds_object = dds_subset_filt
  #contrast_list = list_of_contrasts[[1]]
  
  res_table_unshrunken <- DESeq2::results(dds_object, 
                                          contrast = contrast_list,
                                          parallel = TRUE, 
                                          alpha = padj_cutoff)
  res_table <- DESeq2::lfcShrink(dds_object, 
                                 contrast = contrast_list,
                                 res = res_table_unshrunken, 
                                 type = "ashr"
                                 )
  
  results_data_annot <- as.data.frame(res_table) %>%
    tibble::rownames_to_column("ensembl_id") %>% 
    dplyr::mutate(FoldChange=ifelse(log2FoldChange < 0, -2^(abs(log2FoldChange)), 2^(log2FoldChange))) %>% #adding FoldChange column
    dplyr::left_join(ensemblAnnot, "ensembl_id") %>% # adding annotationl entrez_ids, gene symbols,...
    # dplyr::left_join(normalized_counts_AddedMean, "ensembl_id") %>% # adding normalized counts and average counts per group
    dplyr::arrange(padj) %>% # ordering based on p.adj
    dplyr::select(ensembl_id, contains("_symbol"), contains("description"), entrezgene,
                  log2FoldChange, lfcSE, FoldChange, contains("stat"), pvalue, padj,
                  gene_biotype, chromosome_name,  start_position, end_position, strand) 
  
  results_data_annot_signif <- results_data_annot %>%
    dplyr::filter((!is.na(padj) & (padj < padj_cutoff)) & abs(log2FoldChange) > log2FC_cutoff)
  
  temp_results_summary_df <- data.frame(test = coeff_name,
                                        design = paste(as.character(design(dds_object)), collapse=""),
                                        signif_genes = nrow(results_data_annot_signif),
                                        signif_genes_UP = sum(results_data_annot_signif$log2FoldChange > log2FC_cutoff),
                                        signif_genes_DOWN = sum(results_data_annot_signif$log2FoldChange < log2FC_cutoff),
                                        cutoffs = paste0("(!is.na(padj) & (padj < ", padj_cutoff,")) & abs(log2FoldChange) > ", log2FC_cutoff))
  
  results_list <- list(results_signif = results_data_annot_signif,
                       de_details = temp_results_summary_df,
                       results_all = results_data_annot)
  return(results_list)
}
