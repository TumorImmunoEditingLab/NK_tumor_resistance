######################################
library(ggbreak)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(extrafont)
library(ggrastr)

import::from(.from = hypeR, ggempty)
import::from(.from = here::here("utils/filterDatasets.R"), "filterDatasets", .character_only=TRUE) # used for filtering
import::from(.from = here::here("utils/generateEnsemblAnnotation.R"), "generateEnsemblAnnotation", .character_only=TRUE) # used for filtering
import::from(.from = here::here("utils/generateResults.R"), "meanExprsPerGroup", "generateResults_upd", "generateResults_upd_for_review", .character_only=TRUE) # used for filtering
import::from(.from = here::here("utils/plotDegResults.R"), "plotVolcano", "plotVolcano_repel", .character_only=TRUE)
import::from(.from = GenomicFeatures, makeTxDbFromGFF)
import::from(.from = AnnotationDbi, annot_db_keys = keys, annot_db_select = select)
import::from(.from = DESeq2, .all=TRUE)
import::from(.from = tximport, tximport)
######################################

# load custom functions
mapfun <- function(mousegenes) {
  gns <- mapIds(org.Mm.eg.db, mousegenes, "ENTREZID", "SYMBOL")
  mapped <- select(Orthology.eg.db, gns, "Homo_sapiens","Mus_musculus")
  naind <- is.na(mapped$Homo_sapiens)
  hsymb <- mapIds(org.Hs.eg.db, as.character(mapped$Homo_sapiens[!naind]), "SYMBOL", "ENTREZID")
  out <- data.frame(Mouse_symbol = mousegenes, mapped, Human_symbol = NA)
  out$Human_symbol[!naind] <- hsymb
  out
}

hyp_dots <- function(hyp_obj,
                     top = 20,
                     abrv = 50,
                     size_by = c("genesets", "significance", "none"),
                     pval = 1, 
                     fdr = 1,
                     val = c("fdr", "pval"), 
                     title = "",
                     merge = FALSE) {
  
  stopifnot(is(hyp_obj, "hyp") | is(hyp_obj, "multihyp"))
  
  # Default arguments
  val <- match.arg(val)
  size_by <- match.arg(size_by)
  
  # Handling of multiple signatures
  if (is(hyp_obj, "multihyp")) {
    multihyp_obj <- hyp_obj
    
    # Merge multiple signatures into a single plot
    if (merge) {
      .dots_multi_plot(multihyp_obj$data, top, abrv, size_by, pval, fdr, val, title)
    } 
    # Return a list of plots for each signature
    else {
      mapply(function(hyp_obj, title) {
        
        hyp_dots(hyp_obj,
                 top=top,
                 abrv=abrv,
                 size_by=size_by,
                 pval=pval,
                 fdr=fdr,
                 val=val,
                 title=title)
        
      }, multihyp_obj$data, names(multihyp_obj$data), USE.NAMES=TRUE, SIMPLIFY=FALSE)
    }
  } 
  else {
    .dots_plot(hyp_obj$data, top, abrv, size_by, pval, fdr, val, title)
  }
}

# from "montilab/hypeR" 
.reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv, 
                    scales::log_breaks(base=base), 
                    domain=c(1e-100, Inf))
}

# from "montilab/hypeR" 
.dots_plot <- function(hyp_df,
                       top = 20,
                       abrv = 50,
                       size_by = c("genesets", "significance", "none"),
                       pval_cutoff = 1, 
                       fdr_cutoff = 1,
                       val = c("fdr", "pval"),
                       title = "") {
  
  # Default arguments
  val <- match.arg(val)
  size_by <- match.arg(size_by)
  
  # Subset results
  df <- hyp_df %>%
    dplyr::filter(pval <= pval_cutoff) %>%
    dplyr::filter(fdr <= fdr_cutoff) %>%
    purrr::when(!is.null(top) ~ head(., top), ~ .)
  
  # Handle empty dataframes
  if (nrow(df) == 0) return(ggempty())
  
  # Plotting variables
  df$significance <- df[,val]
  df$size <- 1
  
  if (size_by == "significance") {
    df$size <- df$significance
  }
  if (size_by == "genesets") {
    df$size <- df$geneset
  }
  
  # Order by significance value
  df <- df[order(-df[,val]),]
  
  # Abbreviate labels
  label.abrv <- substr(df$label, 1, abrv)
  if (any(duplicated(label.abrv))) {
    stop("Non-unique labels after abbreviating")
  } else {
    df$label.abrv <- factor(label.abrv, levels=label.abrv)   
  }
  
  if (val == "pval") {
    color.label <- "P-Value"
  }
  if (val == "fdr") {
    color.label <- "FDR"
  }
  
  p <- ggplot(df, aes(x=label.abrv, y=significance, color=significance, size=size)) +
    geom_point() +
    labs(title=title, y=color.label, color=color.label) +
    scale_color_continuous(low="#E53935", high="#114357", guide=guide_colorbar(reverse=TRUE)) +
    coord_flip() +
    scale_y_continuous(trans=.reverselog_trans(10)) +
    geom_hline(yintercept=0.05, linetype="dotted") +
    theme(plot.title=element_text(hjust=0.5),
          axis.title.y=element_blank()) +
    guides(size = guide_legend(order=2))
  
  if (size_by == "none") {
    p <- p + guides(size="none")
  }
  if (size_by == "significance") {
    p <- p + scale_size_continuous(trans=.reverselog_trans(10)) + labs(size="Significance")
  }
  if (size_by == "genesets") {
    p <- p + scale_size_continuous(trans=scales::log10_trans()) + labs(size="Genesets\nSize")
  }
  
  return(p)
}
import::from(.from = here::here("utils/plotDegResults.R"), "plotVolcano", "plotVolcano_repel", .character_only=TRUE)
library(hypeR)

######
# declare options
experiment_name="nk_tum_immunoedit_new_data" # this can be a subset analysis - e.g. just batch1,...
experiment_design="1"  # used only to construct initial dds object
abs_filt_samples=3

# parameters for annotation
biomart_host = "http://nov2020.archive.ensembl.org"
#biomart_host = "https://www.ensembl.org"
biomart_Ens_version = "Ensembl Genes 102"
biomart_dataset="mmusculus_gene_ensembl"


##############################
# Plotting GRCm38 results from MS

# =================================================================================================
# declare options 
options(ggrepel.max.overlaps = Inf)
padj_cutoff <- 0.05
log2FC_cutoff <- 0.58

# actual processign
path_data_folder <- "~/workspace/datasets/Plots to edit/" 
path_data_files <- list.files(path = path_data_folder,
                              pattern = "*differential.tsv",
                              full.names = TRUE)
head_number <- 10


pdf(file = "~/workspace/results/review_additions/volcanos_10.pdf")
for(path_file in path_data_files){
  path_file = "/home/rstudio/workspace/datasets/Plots to edit//balanced_woo_contrast_treatment_WTNK_versus_untreated_for_A_TP1_differential.tsv"
  path_file = "/home/rstudio/workspace/datasets/Plots to edit//balanced_woo_contrast_treatment_WTNK_versus_untreated_for_A_TP2_differential.tsv"
  #path_file = "/home/rstudio/workspace/datasets/Plots to edit//balanced_woo_contrast_treatment_IfngKONK_versus_untreated_for_A_TP1_differential.tsv"
  #path_file = "/home/rstudio/workspace/datasets/Plots to edit//balanced_woo_contrast_treatment_IfngKONK_versus_untreated_for_TP1_D_differential.tsv"
  #path_file = "/home/rstudio/workspace/datasets/Plots to edit//balanced_woo_contrast_treatment_IfngKONK_versus_untreated_for_TP2_A_differential.tsv"
  #path_file = "/home/rstudio/workspace/datasets/Plots to edit//balanced_woo_contrast_treatment_IfngKONK_versus_untreated_for_TP2_D_differential.tsv"
  #path_file = "/home/rstudio/workspace/datasets/Plots to edit//balanced_woo_contrast_treatment_IfngKONK_versus_untreated_for_TP1_D_differential.tsv"
  
  print(path_file)
  data_tmp <- read.table(file = path_file, header = TRUE, sep = "\t")
  data_tmp <- data_tmp %>% 
    dplyr::filter(!is.na(padj))
  data_tmp$signif_DE <- "NO"
  data_tmp$signif_DE[data_tmp$log2FoldChange > 
                       log2FC_cutoff & data_tmp$padj < padj_cutoff] <- "UP"
  data_tmp$signif_DE[data_tmp$log2FoldChange < 
                       -log2FC_cutoff & data_tmp$padj < padj_cutoff] <- "DOWN"
  
  genes_to_label <- c("Ly6a", "Plaat3", "Lgals3bp", "Stat1", "Ccl5")
  
  tmp_genes_to_label <- data_tmp %>% 
    filter(signif_DE == "DOWN") %>%
    arrange(padj) %>%
    head(head_number) %>%
    pull(gene_name)
  genes_to_label <- c(genes_to_label, tmp_genes_to_label)
  
  tmp_genes_to_label <- data_tmp %>% 
    filter(signif_DE == "UP") %>%
    arrange(padj) %>%
    head(head_number) %>%
    pull(gene_name)
  genes_to_label <- c(genes_to_label, tmp_genes_to_label)
  
  tmp_genes_to_label <- data_tmp %>% 
    filter(padj <= padj_cutoff & abs(log2FoldChange) >= log2FC_cutoff) %>%
    arrange(log2FoldChange) %>%
    head(head_number) %>%
    pull(gene_name)
  genes_to_label <- c(genes_to_label, tmp_genes_to_label)
  
  tmp_genes_to_label <- data_tmp %>% 
    filter(padj <= padj_cutoff & abs(log2FoldChange) >= log2FC_cutoff) %>%
    arrange(log2FoldChange) %>%
    tail(head_number) %>%
    pull(gene_name)
  genes_to_label <- c(genes_to_label, tmp_genes_to_label)
  
  genes_to_label <- unique(genes_to_label)
  
  logical_filter <- genes_to_label %in% (data_tmp %>% filter(signif_DE %in% c("DOWN", "UP")) %>% pull(gene_name))
  genes_to_label <- genes_to_label[logical_filter]
  data_tmp$mgi_symbol <- data_tmp$gene_name
  
  titel <- stringr::str_extract(path_file, pattern = "//.*")
  
  if(titel == "//balanced_woo_contrast_treatment_WTNK_versus_IfngKONK_for_D_TP1_differential.tsv"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      scale_y_break(c(20, 21), scales = 0.33)
  }else if(titel == "//balanced_woo_contrast_treatment_WTNK_versus_untreated_for_D_TP1_differential.tsv"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      scale_y_break(c(20, 21), scales = 0.33)
  }else if (titel == "//balanced_woo_contrast_treatment_IfngKONK_versus_untreated_for_A_TP1_differential.tsv"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      ylim(-0.1, 37) +
      scale_x_continuous(breaks = c(-4, 0, 4, 8, 12), limits = c(-6, 12.5))
  }else if (titel == "//balanced_woo_contrast_treatment_IfngKONK_versus_untreated_for_D_TP1_differential.tsv"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      ylim(-0.1, 20) +
      scale_x_continuous(breaks = c(-2.5, 0, 2.5, 5), limits = c(-2.7, 6.4))
  }else if (titel == "//balanced_woo_contrast_treatment_IfngKONK_versus_untreated_for_A_TP2_differential.tsv"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      ylim(-0.1, 33) +
      scale_x_continuous(breaks = c(-6, -3, 0, 3, 6), limits = c(-6.4, 6.1))
  }else if (titel == "//balanced_woo_contrast_treatment_IfngKONK_versus_untreated_for_D_TP2_differential.tsv"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      ylim(-0.1, 16) +
      scale_x_continuous(breaks = c(-5, 0, 5, 10, 15), limits = c(-8, 16))
  }else{
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8))
  }
  
  print(p)
  
  ggsave(filename = paste0("~/workspace/results/review_additions", titel, ".eps"),
         plot = p,
         device="eps"
  )
}
dev.off()


# Enrichments
pdf(file = "~/workspace/results/review_additions/GOBP_Enrichments.pdf")
path_hyp_export <- "~/workspace/results/review_additions/GOBP_Enrichments_tables/"

# HypeR GSEA plots

#gs_hallmark <- hypeR::msigdb_gsets(species = "Mus musculus", category = c("H"), clean=TRUE) 
#gs_C2_kegg <- hypeR::msigdb_gsets(species = "Mus musculus", category = c("C2"), subcategory = "CP:KEGG", clean=TRUE) 
gs_C5_GOBP <- hypeR::msigdb_gsets(species = "Mus musculus", category = c("C5"), subcategory = "GO:BP", clean=TRUE) 
#gs_C5_GOCC <- hypeR::msigdb_gsets(species = "Mus musculus", category = c("C5"), subcategory = "GO:CC", clean=TRUE) 
#gs_C5_GOMF <- hypeR::msigdb_gsets(species = "Mus musculus", category = c("C5"), subcategory = "GO:MF", clean=TRUE) 

for(path_to_deg_file in path_data_files[1:12]){
  #path_to_deg_file <- "/home/rstudio/workspace/datasets/Plots to edit//balanced_woo_contrast_treatment_IfngKONK_versus_untreated_for_A_TP1_differential.tsv"
  #path_to_deg_file <- "/home/rstudio/workspace/datasets/Plots to edit//balanced_woo_contrast_treatment_IfngKONK_versus_untreated_for_D_TP2_differential.tsv"
  
  print(paste0("Processing ", path_to_deg_file))
  deg_table_total <- read.table(file = path_to_deg_file,
                              header = TRUE,
                              sep = "\t")
  
  deg_table_total <- deg_table_total %>% filter(padj < padj_cutoff, abs(log2FoldChange) >= 0.58)
  deg_table_up_reg <- deg_table_total %>% filter(log2FoldChange > 0)
  deg_table_down_reg <- deg_table_total %>% filter(log2FoldChange < 0)
 
  deg_tables_list <- list(
    TOTAL = deg_table_total,
    UP_REG = deg_table_up_reg,
    DOWN_REG = deg_table_down_reg
  )
  
  if(nrow(deg_table_total) == 0){
    print("No differentially expressed genes, skipping this one")
    next
  }
  
 for(deg_table_name in names(deg_tables_list)) {
   print(paste0("Sub processing ", deg_table_name))
   #deg_table_name <- "TOTAL"
   deg_genes <- deg_tables_list[[deg_table_name]]$gene_name
   
   title <- stringr::str_replace(path_to_deg_file, pattern = ".*//(.*)", replacement = "\\1")
   title <- paste0("GOBP_", deg_table_name, "_", title)

   enrichment_C5_GOBP <- hypeR::hypeR(signature = deg_genes, 
                                      genesets = gs_C5_GOBP, 
                                      test = "hypergeometric")
   
   if(!any(enrichment_C5_GOBP$data$fdr < 0.05)){
     print(paste0("No terms passing FDR filter"))
     next
   }
   
   hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOBP,
                       file_path = paste0(path_hyp_export, title, ".xlsx"))
   
   enrichment_C5_GOBP_plot <- hyp_dots(enrichment_C5_GOBP, 
                                       merge = TRUE, 
                                       fdr = 0.05, 
                                       top = 20, 
                                       abrv = 70, 
                                       val = "fdr", 
                                       title = title) +
     theme_bw()
   enrichment_C5_GOBP_plot
   ggsave(filename = paste0(path_hyp_export, title, ".eps"),
          plot = enrichment_C5_GOBP_plot,
          device = "eps"
   )
   
 } 
  
}
dev.off()





pdf(file = "~/workspace/results/review_additions/Ly6AE_enrichments.pdf")
path_hyp_export <- "~/workspace/results/review_additions/Ly6AE_enrichments_tables/"

path_Ly6a_KO_vs_WT <- "/home/rstudio/workspace/datasets/Plots to edit//ly6a_contrast_group_Ly6a_KO_vs_Ly6a_WT_against_intercept_differential.tsv"
Ly6a_KO_vs_WT <- read.table(file = path_Ly6a_KO_vs_WT,
                            header = TRUE,
                            sep = "\t")
Ly6a_KO_vs_WT <- Ly6a_KO_vs_WT %>% filter(padj <= padj_cutoff,
                                          abs(log2FoldChange) >= 0.58)
enrichment_C5_GOBP <- hypeR::hypeR(signature = Ly6a_KO_vs_WT$gene_name, 
                                      genesets = gs_C5_GOBP, 
                                      test = "hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOBP,
                    file_path = paste0(path_hyp_export, "GOBP_Ly6a_KO_vs_WT.xlsx"))
enrichment_C5_GOBP_plot <- hyp_dots(enrichment_C5_GOBP, 
                                       merge = TRUE, 
                                       fdr = 0.05, 
                                       top = 20, 
                                       abrv = 70, 
                                       val = "fdr", 
                                       title = "GOBP Ly6a_KO_vs_WT") +
  theme_bw()
enrichment_C5_GOBP_plot
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOBP_Ly6a_KO_vs_WT.eps"),
       plot = enrichment_C5_GOBP_plot,
       device="eps"
)


enrichment_C5_GOCC <- hypeR::hypeR(signature = Ly6a_KO_vs_WT$gene_name, 
                                      genesets = gs_C5_GOCC, 
                                      test = "hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOCC,
                    file_path = paste0(path_hyp_export, "GOCC_Ly6a_KO_vs_WT.xlsx"))
enrichment_C5_GOCC_plot <- hyp_dots(enrichment_C5_GOCC, 
                                       merge = TRUE, 
                                       fdr = 0.05, 
                                       top = 20, 
                                       abrv = 70, 
                                       val = "fdr", 
                                       title = "GOCC Ly6a_KO_vs_WT") +
  theme_bw()
enrichment_C5_GOCC_plot
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOCC_Ly6a_KO_vs_WT.eps"),
       plot = enrichment_C5_GOCC_plot,
       device="eps"
)





# 
# Ly6a_KO_vs_WT_up <- Ly6a_KO_vs_WT %>% filter(log2FoldChange > 0)
# Ly6a_KO_vs_WT_down <- Ly6a_KO_vs_WT %>% filter(log2FoldChange < 0)
# 
# enrichment_C5_GOBP_up <- hypeR::hypeR(signature = Ly6a_KO_vs_WT_up$gene_name, 
#                                       genesets = gs_C5_GOBP, 
#                                       test = "hypergeometric")
# hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOBP_up,
#                     file_path = paste0(path_hyp_export, "GOBP_Ly6a_KO_vs_WT_up.xlsx"))
# enrichment_C5_GOBP_plot_up <- hyp_dots(enrichment_C5_GOBP_up, 
#                                        merge = TRUE, 
#                                        fdr = 0.05, 
#                                        top = 20, 
#                                        abrv = 70, 
#                                        val = "fdr", 
#                                        title = "GOBP Ly6a_KO_vs_WT_up") +
#   theme_bw()
# enrichment_C5_GOBP_plot_up
# ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOBP_Ly6a_KO_vs_WT_up.eps"),
#        plot = enrichment_C5_GOBP_plot_up,
#        device="eps"
# )
# 
# enrichment_C5_GOBP_down <- hypeR::hypeR(signature = Ly6a_KO_vs_WT_down$gene_name, 
#                                         genesets = gs_C5_GOBP, 
#                                         test = "hypergeometric")
# hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOBP_down,
#                     file_path = paste0(path_hyp_export, "GOBP_Ly6a_KO_vs_WT_down.xlsx"))
# enrichment_C5_GOBP_plot_down <- hyp_dots(enrichment_C5_GOBP_down, 
#                                          merge = TRUE, 
#                                          fdr = 0.05, 
#                                          top = 20, 
#                                          abrv = 70, 
#                                          val = "fdr", 
#                                          title = "GOBP Ly6a_KO_vs_WT_down") +
#   theme_bw()
# enrichment_C5_GOBP_plot_down
# ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOBP_Ly6a_KO_vs_WT_down.eps"),
#        plot = enrichment_C5_GOBP_plot_down,
#        device="eps"
# )
# 
# enrichment_C5_GOCC_up <- hypeR::hypeR(signature = Ly6a_KO_vs_WT_up$gene_name, 
#                                       genesets = gs_C5_GOCC, 
#                                       test = "hypergeometric")
# hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOCC_up,
#                     file_path = paste0(path_hyp_export, "GOCC_Ly6a_KO_vs_WT_up.xlsx"))
# enrichment_C5_GOCC_plot_up <- hyp_dots(enrichment_C5_GOCC_up, 
#                                        merge = TRUE, 
#                                        fdr = 0.05, 
#                                        top = 20, 
#                                        abrv = 70, 
#                                        val = "fdr", 
#                                        title = "GOCC Ly6a_KO_vs_WT_up") +
#   theme_bw()
# enrichment_C5_GOCC_plot_up
# ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOCC_Ly6a_KO_vs_WT_up.eps"),
#        plot = enrichment_C5_GOCC_plot_up,
#        device="eps"
# )
# 
# enrichment_C5_GOCC_down <- hypeR::hypeR(signature = Ly6a_KO_vs_WT_down$gene_name, 
#                                         genesets = gs_C5_GOCC, 
#                                         test = "hypergeometric")
# hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOCC_down,
#                     file_path = paste0(path_hyp_export, "GOCC_Ly6a_KO_vs_WT_down.xlsx"))
# enrichment_C5_GOCC_plot_down <- hyp_dots(enrichment_C5_GOCC_down, 
#                                          merge = TRUE, 
#                                          fdr = 0.05, 
#                                          top = 20, 
#                                          abrv = 70, 
#                                          val = "fdr", 
#                                          title = "GOCC Ly6a_KO_vs_WT_down") +
#   theme_bw()
# enrichment_C5_GOCC_plot_down
# ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOCC_Ly6a_KO_vs_WT_down.eps"),
#        plot = enrichment_C5_GOCC_plot_down,
#        device="eps"
# )


#gs_hallmark_h <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("H"), clean=TRUE) 
#gs_C2_kegg_h <-  hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:KEGG", clean=TRUE) 
gs_C5_GOBP_h <-  hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:BP", clean=TRUE) 
gs_C5_GOCC_h <-  hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:CC", clean=TRUE) 
#gs_C5_GOMF_h <-  hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:MF", clean=TRUE) 

path_Ly6e_KO_vs_WT <- "/home/rstudio/workspace/datasets/Plots to edit//human_contrast_group_LY6E_KO_vs_LY6E_WT_against_intercept_differential.tsv"
Ly6e_KO_vs_WT <- read.table(file = path_Ly6e_KO_vs_WT,
                            header = TRUE,
                            sep = "\t")
Ly6e_KO_vs_WT <- Ly6e_KO_vs_WT %>% filter(padj <= padj_cutoff, abs(log2FoldChange) >= 0.58)

enrichment_C5_GOBP <- hypeR::hypeR(signature = Ly6e_KO_vs_WT$gene_name, 
                                      genesets = gs_C5_GOBP_h, 
                                      test = "hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOBP,
                    file_path = paste0(path_hyp_export, "GOBP_Ly6e_KO_vs_WT.xlsx"))
enrichment_C5_GOBP_plot <- hyp_dots(enrichment_C5_GOBP, 
                                       merge = TRUE, 
                                       fdr = 0.05, 
                                       top = 20, 
                                       abrv = 70, 
                                       val = "fdr", 
                                       title = "GOBP Ly6e_KO_vs_WT") +
  theme_bw()
enrichment_C5_GOBP_plot
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOBP_Ly6e_KO_vs_WT.eps"),
       plot = enrichment_C5_GOBP_plot,
       device="eps"
)

enrichment_C5_GOCC <- hypeR::hypeR(signature = Ly6e_KO_vs_WT$gene_name, 
                                      genesets = gs_C5_GOCC_h, 
                                      test = "hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOCC,
                    file_path = paste0(path_hyp_export, "GOCC_Ly6e_KO_vs_WT.xlsx"))
enrichment_C5_GOCC_plot <- hyp_dots(enrichment_C5_GOCC, 
                                       merge = TRUE, 
                                       fdr = 0.05, 
                                       top = 20, 
                                       abrv = 70, 
                                       val = "fdr", 
                                       title = "GOCC Ly6e_KO_vs_WT") +
  theme_bw()
enrichment_C5_GOCC_plot
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOCC_Ly6e_KO_vs_WT.eps"),
       plot = enrichment_C5_GOCC_plot,
       device="eps"
)



Ly6e_KO_vs_WT_up <- Ly6e_KO_vs_WT %>% filter(log2FoldChange > 0)
Ly6e_KO_vs_WT_down <- Ly6e_KO_vs_WT %>% filter(log2FoldChange < 0)

enrichment_C5_GOBP_up <- hypeR::hypeR(signature = Ly6e_KO_vs_WT_up$gene_name, 
                                      genesets = gs_C5_GOBP_h, 
                                      test = "hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOBP_up,
                    file_path = paste0(path_hyp_export, "GOBP_Ly6e_KO_vs_WT_up.xlsx"))
enrichment_C5_GOBP_plot_up <- hyp_dots(enrichment_C5_GOBP_up, 
                                       merge = TRUE, 
                                       fdr = 0.05, 
                                       top = 20, 
                                       abrv = 70, 
                                       val = "fdr", 
                                       title = "GOBP Ly6e_KO_vs_WT_up") +
  theme_bw()
enrichment_C5_GOBP_plot_up
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOBP_Ly6e_KO_vs_WT_up.eps"),
       plot = enrichment_C5_GOBP_plot_up,
       device="eps"
)

enrichment_C5_GOBP_down <- hypeR::hypeR(signature = Ly6e_KO_vs_WT_down$gene_name, 
                                        genesets = gs_C5_GOBP_h, 
                                        test = "hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOBP_down,
                    file_path = paste0(path_hyp_export, "GOBP_Ly6e_KO_vs_WT_down.xlsx"))
enrichment_C5_GOBP_plot_down <- hyp_dots(enrichment_C5_GOBP_down, 
                                         merge = TRUE, 
                                         fdr = 0.05, 
                                         top = 20, 
                                         abrv = 70, 
                                         val = "fdr", 
                                         title = "GOBP Ly6e_KO_vs_WT_down") +
  theme_bw()
enrichment_C5_GOBP_plot_down
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOBP_Ly6e_KO_vs_WT_down.eps"),
       plot = enrichment_C5_GOBP_plot_down,
       device="eps"
)


enrichment_C5_GOCC_up <- hypeR::hypeR(signature = Ly6e_KO_vs_WT_up$gene_name, 
                                      genesets = gs_C5_GOCC_h, 
                                      test = "hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOCC_up,
                    file_path = paste0(path_hyp_export, "GOCC_Ly6e_KO_vs_WT_up.xlsx"))
enrichment_C5_GOCC_plot_up <- hyp_dots(enrichment_C5_GOCC_up, 
                                       merge = TRUE, 
                                       fdr = 0.05, 
                                       top = 20, 
                                       abrv = 70, 
                                       val = "fdr", 
                                       title = "GOCC Ly6e_KO_vs_WT_up") +
  theme_bw()
enrichment_C5_GOCC_plot_up
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOCC_Ly6e_KO_vs_WT_up.eps"),
       plot = enrichment_C5_GOCC_plot_up,
       device="eps"
)


enrichment_C5_GOCC_down <- hypeR::hypeR(signature = Ly6e_KO_vs_WT_down$gene_name, 
                                        genesets = gs_C5_GOCC_h, 
                                        test = "hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOCC_down,
                    file_path = paste0(path_hyp_export, "GOCC_Ly6e_KO_vs_WT_down.xlsx"))
enrichment_C5_GOCC_plot_down <- hyp_dots(enrichment_C5_GOCC_down, 
                                         merge = TRUE, 
                                         fdr = 0.05, 
                                         top = 20, 
                                         abrv = 70, 
                                         val = "fdr", 
                                         title = "GOCC Ly6e_KO_vs_WT_down") +
  theme_bw()
enrichment_C5_GOCC_plot_down
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOCC_Ly6e_KO_vs_WT_down.eps"),
       plot = enrichment_C5_GOCC_plot_down,
       device="eps"
)

# check what are the orthologs of mice genes in human and do GOBP
library(Orthology.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

Human_orthologs_ly6a <- mapfun(Ly6a_KO_vs_WT$gene_name)
Human_orthologs_ly6a <- Human_orthologs_ly6a$Human_symbol
Human_orthologs_ly6a <- Human_orthologs_ly6a[!is.na(Human_orthologs_ly6a)]

enrichment_C5_GOBP <- hypeR::hypeR(signature = Human_orthologs_ly6a,
                                      genesets = gs_C5_GOBP_h,
                                      test="hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOBP,
                    file_path = paste0(path_hyp_export, "GOBP_Ly6a_KO_vs_WT_total-regulated_HUMAN_ORTHOLOGUES.xlsx"))
enrichment_C5_GOBP_plot <- hyp_dots(enrichment_C5_GOBP,
                                       merge=TRUE,
                                       fdr=0.05,
                                       top = 20,
                                       abrv=70,
                                       val="fdr",
                                       title="GOBP Ly6a_KO_vs_WT_total-regulated HUMAN ORTHOLOGUES") +
  theme_bw()
enrichment_C5_GOBP_plot
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOBP_Ly6a_KO_vs_WT_total-regulated_HUMAN_ORTHOLOGUES.eps"),
       plot = enrichment_C5_GOBP_plot,
       device="eps"
)

dev.off()


###### Heat map for tumore_only TP1 + WTNK TP1 + WTNK TP2
padj_cutoff <- 0.05
log2FC_cutoff <- 0.58

rnaseq_deseq <- readRDS("~/workspace/datasets/Plots to edit/rnaseq_deseq_balanced_woo_deseq_data_set.rds")

# 
# rnaseq_deseq$group <- droplevels(rnaseq_deseq$group)
# design(rnaseq_deseq) <- as.formula("~ group")
vsd_subset_filt <- DESeq2::vst(rnaseq_deseq, blind = FALSE) 
vsd_subset_filt <- vsd_subset_filt[, vsd_subset_filt$group %in% c("TP1_Tumor_A", "TP1_Tumor_A_WTNK", "TP2_Tumor_A_WTNK")]

path_file = "/home/rstudio/workspace/datasets/Plots to edit//balanced_woo_contrast_treatment_WTNK_versus_untreated_for_A_TP1_differential.tsv"
data_tmp <- read.table(file = path_file, header = TRUE, sep = "\t")
DE_gene_names <- data_tmp %>% dplyr::filter(!is.na(padj)) %>%
  filter(abs(log2FoldChange) >= log2FC_cutoff & padj <= padj_cutoff) %>%
  pull(gene_id)

path_file = "/home/rstudio/workspace/datasets/Plots to edit//balanced_woo_contrast_treatment_WTNK_versus_untreated_for_A_TP2_differential.tsv"
data_tmp <- read.table(file = path_file, header = TRUE, sep = "\t")
DE_gene_names <- c(DE_gene_names,
                   data_tmp %>% dplyr::filter(!is.na(padj)) %>%
                     filter(abs(log2FoldChange) >= log2FC_cutoff & padj <= padj_cutoff) %>%
                     pull(gene_id)
)
DE_gene_names <- unique(DE_gene_names)

# rnaseq_deseq$sample <- droplevels(rnaseq_deseq$sample)
# rnaseq_deseq$cell_line <- droplevels(rnaseq_deseq$cell_line)
# rnaseq_deseq$treatment <- droplevels(rnaseq_deseq$rnaseq_deseq$treatment)
# 
vsd_subset_filt_assay <- SummarizedExperiment::assay(vsd_subset_filt) 
vsd_subset_filt_assay <- vsd_subset_filt_assay[row.names(vsd_subset_filt_assay) %in% DE_gene_names, ]

metadata_heatmap <- as.data.frame(SummarizedExperiment::colData(vsd_subset_filt))

metadata_heatmap <- metadata_heatmap %>%
  dplyr::select(treatment, time) 

color.scheme <- rev(RColorBrewer::brewer.pal(8,"RdBu")) # generate the color scheme to use

ann_colors = list(
  treatment = c( untreated = "#005f73", WTNK = "#0a9396"),
  time = c(TP1 = "#DD3344", TP2 = "#FF9F1C")
)

#row_annot_symbols <- data.frame(genes = c("Plaat3", "Ly6A"))
#row.names(row_annot_symbols) <- c("ENSMUSG00000060675", "ENSMUSG00000075602")

heatmap_corrected <- pheatmap::pheatmap(vsd_subset_filt_assay,
                                        main = "Heatmap of signif. DEG - Batch Corrected",
                                        scale = "row",
                                        annotation_col = metadata_heatmap,
                                        annotation_colors = ann_colors,
 #                                       annotation_row = row_annot_symbols,
                                        show_colnames = FALSE,
                                        show_rownames = FALSE,
                                        cluster_cols = FALSE,
                                        #cluster_rows = counts_deg_ord_row_cor_hclust,
                                        color = color.scheme,
                                        fontsize = 10, fontsize_row = 10 #height=10, cellwidth = 11, cellheight = 11
) 
ggsave("~/workspace/results/review_additions/heatmap_TP1_2_Tumor_A_WTNK_vs_WTNK_A_TP1.eps",
       heatmap_corrected,
       device = "eps", width = 210, height = 300, dpi = 320, units = "mm"
       )

ggsave("~/workspace/results/review_additions/heatmap_TP1_2_Tumor_A_WTNK_vs_WTNK_A_TP1.pdf",
       heatmap_corrected,
       device = "pdf", width = 210, height = 300, dpi = 320, units = "mm"
)






#############################
# load the raw data
rnaseq_metadata_df_raw <- readr::read_csv2(file = "~/workspace/data_for_review/exp_for_review_metadata.csv")  # 133 samples; 42 samples from other project?
quantseq_files <- "~/workspace/data_for_review/salmon_quant/"

rnaseq_metadata_df_raw

metadata$absolute_quant_files_path
metadata$sample_name

metadata <- rnaseq_metadata_df_raw

quant_dirs <- list.dirs("~/workspace/data_for_review/salmon_quant/", recursive = F)
quant_dirs <- paste0(quant_dirs,"/quant.sf")
quant_dirs_df <- data.frame(full_path = quant_dirs, sample_name = str_extract(quant_dirs, pattern = "(X..)"))

metadata$sample_name

metadata <- merge.data.frame(metadata, quant_dirs_df, by = "sample_name")
metadata$timepoint_cell_harvesting <- str_replace_all(metadata$timepoint_cell_harvesting, "timepoint ", "TP")

metadata$sample_name_used <- paste0(metadata$condition, "_", metadata$timepoint_cell_harvesting, "_", metadata$cell_line_label)

annotation_dir <- "~/workspace/datasets/annotations/"
gtf_file <- file.path(annotation_dir, "Mus_musculus.GRCm38.102.gtf") # needed only once to generate tx2gene.tsv
tx2gene_file <- file.path(annotation_dir, "EnsDb.Mmusculus.v102_tx2gene.tsv")
tx2gene_file_md5sum <- file.path(annotation_dir, "EnsDb.Mmusculus.v102_tx2gene.md5sum")

if (!file.exists(tx2gene_file)) {
  message(tx2gene_file, " does not exist: generating!")
  txDB <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format="gtf")
  
  # AnnotationDbi::keys
  tx_name <- annot_db_keys(txDB, keytype="TXNAME")
  
  # AnnotationDbi::select
  tx2gene_gtf <- annot_db_select(txDB,
                                 keys = tx_name,
                                 columns="GENEID",
                                 keytype = "TXNAME")
  
  
  # saving tx2gene into destdir
  write.table(tx2gene_gtf,
              file = tx2gene_file,
              sep = "\t",
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE)
  
  tx2gene_md5sum <- digest::digest(tx2gene_file, file=TRUE, algo="md5")
  write.table(paste(tx2gene_md5sum, basename(tx2gene_file), collapse = "\t"),
              file = tx2gene_file_md5sum,
              sep = "\t",
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE)
  
} else {
  message(tx2gene_file, " exists: loading!")
  tx2gene <- readr::read_tsv(tx2gene_file,
                             col_names = c("tx_id", "gene_id"),
                             show_col_types = FALSE)
}

output_dir <- "~/workspace/RDSs/"
precomputed_objects_filename <- paste0(experiment_name, "_dds_objects.RData")
precomputed_objects_file <- file.path(output_dir, precomputed_objects_filename)

precomputed_transf_objects_filename <- paste0(experiment_name, "_log2_vsd_filt_objects.RData")
precomputed_transf_objects_file <- file.path(output_dir, precomputed_transf_objects_filename)

if(!file.exists(precomputed_objects_file)){
  cat("RNA objects file '", precomputed_objects_filename, "' does not exist in the OUTPUT_DIR...\n")
  cat("generating", precomputed_objects_filename, "\n")
  # Generating dds and filtered dds objects ----
  #require(DESeq2)
  # here there is a single sample so we use ~1.
  # expect a warning that there is only a single sample...
  
  # generating dds object ----
  # txi_object <- generateTxi(metadata_object=metadata,
  #                           column_names=c("absolute_quant_files_path", "filenames"),
  #                           tx2gene = tx2gene)
  
  quant_files <- metadata$full_path
  names(quant_files) <- metadata$sample_name
  rownames(metadata) <- metadata$sample_name
  
  txi_object <- tximport::tximport(quant_files,
                                   type="salmon",
                                   tx2gene = tx2gene,
                                   ignoreTxVersion = TRUE)
  
  sum(rownames(metadata) == colnames(txi_object$counts)) # check if names match between metadata and counts data
  cat("... txi_object\n")
  
  dds <- DESeq2::DESeqDataSetFromTximport(txi_object, 
                                          colData = metadata, 
                                          design = as.formula(paste0("~", experiment_design)))
  
  cat("... dds\n")
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::DESeq(dds, minReplicatesForReplace=Inf) # do not replace outliers based on replicates
  
  # filtering dds
  cat("... dds filtering\n")
  dds_filt <- filterDatasets(dds, abs_filt = TRUE, 
                             abs_filt_samples = abs_filt_samples) # at least in N samples, which is a smallest group size
  dds_filt <- DESeq2::estimateSizeFactors(dds_filt)
  
  # Generate annotation file using biomaRt
  #FIXME:
  # Issues with Curl when using https:// and http:// throws a warning:
  #  Warning: Ensembl will soon enforce the use of https.
  #   Ensure the 'host' argument includes "https://"
  cat("... fetching annotation from biomart\n")
  ensemblAnnot <- generateEnsemblAnnotation(ensembl_ids = rownames(dds),
                                            host=biomart_host,
                                            version=biomart_Ens_version,
                                            dataset=biomart_dataset)
  
  cat("saving...", precomputed_objects_file, "\n")
  cat("...into:", output_dir, "\n")
  save(dds, dds_filt, ensemblAnnot,  file = precomputed_objects_file)
  
} else{
  cat("RNA objects file '", precomputed_objects_filename, "' exist...loading\n")
  load(precomputed_objects_file)
}


### Load GOs
gs_C5_GOBP <- hypeR::msigdb_gsets(species = "Mus musculus", category = c("C5"), subcategory = "GO:BP", clean=TRUE) 
gs_C5_GOCC <- hypeR::msigdb_gsets(species = "Mus musculus", category = c("C5"), subcategory = "GO:CC", clean=TRUE) 
#gs_C5_GOMF <- hypeR::msigdb_gsets(species = "Mus musculus", category = c("C5"), subcategory = "GO:MF", clean=TRUE) 



#### analysis of ly6aKO vs WT
abs_filt_samples=3
padj_cutoff = 0.05
log2FC_cutoff = 0.58 #(FC=1.5); log2FC=1.0 # (FC=2)
var_expl_needed <- 0.6   

if(exists("deg_design")) {rm(deg_design)}
deg_design <- as.formula("~ condition")

dds_subset <- dds_filt[ , dds_filt$original_experiment %in% "EXP 15.4"]  

dds_subset$condition <- factor(dds_subset$condition, levels = c("WT_Ly6a", "KO_Ly6a"))

dds_subset$condition <- droplevels(dds_subset$condition)

table(dds_subset$condition)

dds_subset_filt <- filterDatasets(dds_subset, 
                                  abs_filt = TRUE, 
                                  abs_filt_samples = abs_filt_samples)

design(dds_subset_filt) <- deg_design
dds_subset_filt <- DESeq2::estimateSizeFactors(dds_subset_filt)
dds_subset_filt <- DESeq2::DESeq(dds_subset_filt) # do not replace outliers based on replicates

log2_norm_subset_filt <- DESeq2::normTransform(dds_subset_filt)
vsd_subset_filt <- DESeq2::vst(dds_subset_filt, blind = FALSE) 

resultsNames(dds_subset_filt)
#debug(generateResults_upd)

#results(dds_subset_filt, contrast=c("condition", "KO_Ly6a", "WT_Ly6a"))
coeff_name <- "condition_KO_Ly6a_vs_WT_Ly6a"
deg_condition_KO_Ly6a_vs_WT_Ly6a_results <- generateResults_upd(dds_object = dds_subset_filt, 
                                                         coeff_name = coeff_name,
                                                         cond_numerator = "KO_Ly6a", 
                                                         cond_denominator = "WT_Ly6a",
                                                         cond_variable = "condition")

openxlsx::write.xlsx(deg_condition_KO_Ly6a_vs_WT_Ly6a_results, 
                     file = file.path("~/workspace/results/review_additions", paste0("deg_", coeff_name, ".xlsx")))

deg_condition_KO_Ly6a_vs_WT_Ly6a_results
#genes_of_interest <- c("Ly6a", "Plaat3", "Lgals3bp", "Stat1", "Ccl5")

data_tmp <- deg_condition_KO_Ly6a_vs_WT_Ly6a_results$results_all
data_tmp$signif_DE <- "NO"
data_tmp$signif_DE[data_tmp$log2FoldChange > log2FC_cutoff & 
                     data_tmp$padj < padj_cutoff] <- "UP"
data_tmp$signif_DE[data_tmp$log2FoldChange < -log2FC_cutoff & 
                     data_tmp$padj < padj_cutoff] <- "DOWN"

data_tmp$gene_name <- data_tmp$mgi_symbol

genes_to_label <- c("Ly6a", "Plaat3", "Lgals3bp", "Stat1", "Ccl5")
head_number <- 10

tmp_genes_to_label <- data_tmp %>% 
  filter(signif_DE == "DOWN") %>%
  arrange(padj) %>%
  head(head_number) %>%
  pull(gene_name)
genes_to_label <- c(genes_to_label, tmp_genes_to_label)

tmp_genes_to_label <- data_tmp %>% 
  filter(signif_DE == "UP") %>%
  arrange(padj) %>%
  head(head_number) %>%
  pull(gene_name)
genes_to_label <- c(genes_to_label, tmp_genes_to_label)

tmp_genes_to_label <- data_tmp %>% 
  filter(padj <= padj_cutoff & abs(log2FoldChange) >= log2FC_cutoff) %>%
  arrange(log2FoldChange) %>%
  head(head_number) %>%
  pull(gene_name)
genes_to_label <- c(genes_to_label, tmp_genes_to_label)

tmp_genes_to_label <- data_tmp %>% 
  filter(padj <= padj_cutoff & abs(log2FoldChange) >= log2FC_cutoff) %>%
  arrange(log2FoldChange) %>%
  tail(head_number) %>%
  pull(gene_name)
genes_to_label <- c(genes_to_label, tmp_genes_to_label)

genes_to_label <- unique(genes_to_label)

logical_filter <- genes_to_label %in% (data_tmp %>% filter(signif_DE %in% c("DOWN", "UP")) %>% pull(gene_name))
genes_to_label <- genes_to_label[logical_filter]
data_tmp$mgi_symbol <- data_tmp$gene_name

TpNK_tp2_vs_Tonly_tp1_volcano_plot_1 <- plotVolcano(dds_results_obj = deg_condition_KO_Ly6a_vs_WT_Ly6a_results$results_all, 
                                                    genes_of_interest = genes_to_label, 
                                                    plot_title = "deg_condition_KO_Ly6a_vs_WT_Ly6a_results")
ggsave(TpNK_tp2_vs_Tonly_tp1_volcano_plot_1, filename = paste0("~/workspace/results/review_additions/deg_", coeff_name,"_volcanoPlot.eps"), 
       device = "eps", width = 7, height = 7, units = "in")

enrichment_C5_GOBP_up <- hypeR::hypeR(signature = deg_condition_KO_Ly6a_vs_WT_Ly6a_results$results_signif %>% filter(log2FoldChange > 0) %>% pull(mgi_symbol), 
                                      genesets = gs_C5_GOBP, 
                                      test = "hypergeometric")
# hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOBP_up,
#                     file_path = paste0(path_hyp_export, "GOBP_Ly6a_KO_vs_WT_up.xlsx"))
enrichment_C5_GOBP_plot_up <- hyp_dots(enrichment_C5_GOBP_up, 
                                       merge = TRUE, 
                                       fdr = 0.05, 
                                       top = 20, 
                                       abrv = 70, 
                                       val = "fdr", 
                                       title = "GOBP Ly6a_KO_vs_WT_up") +
  theme_bw()
ggsave(enrichment_C5_GOBP_plot_up, filename = paste0("~/workspace/results/review_additions/deg_", coeff_name,"_GOBP_UP.eps"),
       device = "eps", width = 7, height = 7, units = "in")

enrichment_C5_GOBP_down <- hypeR::hypeR(signature = deg_condition_KO_Ly6a_vs_WT_Ly6a_results$results_signif %>% filter(log2FoldChange < 0) %>% pull(mgi_symbol), 
                                      genesets = gs_C5_GOBP, 
                                      test = "hypergeometric")
# hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOBP_down,
#                     file_path = paste0(path_hyp_export, "GOBP_Ly6a_KO_vs_WT_up.xlsx"))
enrichment_C5_GOBP_plot_down <- hyp_dots(enrichment_C5_GOBP_down, 
                                         merge = TRUE, 
                                         fdr = 0.05, 
                                         top = 20, 
                                         abrv = 70, 
                                         val = "fdr", 
                                         title = "GOBP Ly6a_KO_vs_WT_down") +
  theme_bw()
ggsave(enrichment_C5_GOBP_plot_down, filename = paste0("~/workspace/results/review_additions/deg_", coeff_name,"_GOBP_DOWN.eps"),
       device = "eps", width = 7, height = 7, units = "in")

#################################################################################
abs_filt_samples=3
padj_cutoff = 0.05
log2FC_cutoff = 0.58 #(FC=1.5); log2FC=1.0 # (FC=2)
var_expl_needed <- 0.6   

if(exists("deg_design")) {rm(deg_design)}
deg_design <- as.formula("~ cell_line_label * condition * timepoint_cell_harvesting ")

dds_subset <- dds_filt[ , dds_filt$original_experiment %in% "EXP 9.11"]  
dds_subset <- dds_subset[ , !dds_subset$condition %in% "Tumor+Prf1_KO_NK"] 
dds_subset$condition <- str_replace_all(dds_subset$condition, pattern = "\\+", replacement = "_")
dds_subset <- dds_subset[ , !dds_subset$timepoint_cell_harvesting %in% "TP0"] 

dds_subset$condition <- factor(dds_subset$condition, levels = c("Tumor_only", "Tumor_Ifng_KO_NK", "Tumor_WT_NK"))
table(dds_subset$condition)
dds_subset$cell_line_label <- factor(dds_subset$cell_line_label , levels = c("A", "D"))
table(dds_subset$cell_line_label)
dds_subset$timepoint_cell_harvesting <- factor(dds_subset$timepoint_cell_harvesting , levels = c("TP1", "TP2"))
table(dds_subset$timepoint_cell_harvesting)

dds_subset_filt <- filterDatasets(dds_subset, 
                               abs_filt = TRUE, 
                               abs_filt_samples = abs_filt_samples)
design(dds_subset_filt) <- deg_design

dds_subset_filt <- DESeq2::estimateSizeFactors(dds_subset_filt)
# dds_subset_filt <- DESeq2::DESeq(dds_subset_filt, 
#                                  test = "LRT", 
#                                  reduced = as.formula("~ cell_line_label * timepoint_cell_harvesting ")) 
dds_subset_filt <- DESeq2::DESeq(dds_subset_filt)
#dds_subset_filt <- DESeq2::DESeq(dds_subset_filt, test = c("Wald"))

# resultsNames(dds_subset_filt)
# results(dds_subset_filt, contrast = list(c("condition_Tumor_Ifng_KO_NK_vs_Tumor_only"), c("Intercept")))
# results(dds_subset_filt, contrast = list(c("condition_Tumor_Ifng_KO_NK_vs_Tumor_only")))

list_of_contrasts <- list(
  IfngKONK_versus_untreated_for_A_TP1 = list(c("condition_Tumor_Ifng_KO_NK_vs_Tumor_only")),
  IfngKONK_versus_untreated_for_D_TP1 = list(c("condition_Tumor_Ifng_KO_NK_vs_Tumor_only", "cell_line_labelD.conditionTumor_Ifng_KO_NK")),
  IfngKONK_versus_untreated_for_A_TP2 = list(c("condition_Tumor_Ifng_KO_NK_vs_Tumor_only", "conditionTumor_Ifng_KO_NK.timepoint_cell_harvestingTP2")),
  IfngKONK_versus_untreated_for_D_TP2 = list(c("condition_Tumor_Ifng_KO_NK_vs_Tumor_only", "cell_line_labelD.conditionTumor_Ifng_KO_NK.timepoint_cell_harvestingTP2")),
  WTNK_versus_untreated_for_A_TP1 = list(c("condition_Tumor_WT_NK_vs_Tumor_only")),
  WTNK_versus_untreated_for_D_TP1 = list(c("condition_Tumor_WT_NK_vs_Tumor_only", "cell_line_labelD.conditionTumor_WT_NK")),
  WTNK_versus_untreated_for_A_TP2 = list(c("condition_Tumor_WT_NK_vs_Tumor_only", "conditionTumor_WT_NK.timepoint_cell_harvestingTP2")),
  WTNK_versus_untreated_for_D_TP2 = list(c("condition_Tumor_WT_NK_vs_Tumor_only", "cell_line_labelD.conditionTumor_WT_NK.timepoint_cell_harvestingTP2"))
)


# cell_line_D_vs_A,cell_lineD.treatmentIfngKONK.timeTP2	intercept

# coeff_name <- resultsNames(dds_subset_filt)[1]
# dds_tmp_filt_results <- generateResults_upd(dds_object = dds_subset_filt, 
#                                             coeff_name = coeff_name,
#                                             cond_numerator = comparison_name[1], 
#                                             cond_denominator = comparison_name[2],
#                                             cond_variable = "sample_name_used")

for(contrast_list_name in names(list_of_contrasts)){
  # contrast_list_name <- names(list_of_contrasts)[1]
  contrast_list <- list_of_contrasts[[contrast_list_name]]
  deg_results <- generateResults_upd_for_review(dds_object = dds_subset_filt,
                                                coeff_name = contrast_list_name,
                                                contrast_list = contrast_list, 
                                                padj_cutoff = padj_cutoff,
                                                log2FC_cutoff = log2FC_cutoff) 
  
  openxlsx::write.xlsx(deg_results, 
                       file = file.path("~/workspace/results/review_additions", paste0("deg_", contrast_list_name, ".xlsx")))
  
  
}


path_data_folder <- "~/workspace/results/review_additions/"
path_data_files <- list.files(path = path_data_folder,
                              pattern = "deg_(Ifng)|(WTNK).*.xlsx",
                              full.names = TRUE)
head_number <- 10

pdf(file = "~/workspace/results/review_additions/volcanos_10.pdf")
for(path_file in path_data_files){
  #path_file <- "/home/rstudio/workspace/results/review_additions//deg_sample_name_used_Tumor_Ifng_KO_NK_TP1_A_vs_Tumor_only_TP1_A.xlsx"
  
  print(path_file)
  data_tmp <- openxlsx::read.xlsx(xlsxFile = path_file,sheet = 3)
  data_tmp <- data_tmp %>% 
    dplyr::filter(!is.na(padj))
  data_tmp$signif_DE <- "NO"
  data_tmp$signif_DE[data_tmp$log2FoldChange > log2FC_cutoff & data_tmp$padj < padj_cutoff] <- "UP"
  data_tmp$signif_DE[data_tmp$log2FoldChange < -log2FC_cutoff & data_tmp$padj < padj_cutoff] <- "DOWN"
  
  genes_to_label <- c("Ly6a", "Plaat3", "Lgals3bp", "Stat1", "Ccl5")
  
  tmp_genes_to_label <- data_tmp %>% 
    filter(signif_DE == "DOWN") %>%
    arrange(padj) %>%
    head(head_number) %>%
    pull(mgi_symbol)
  genes_to_label <- c(genes_to_label, tmp_genes_to_label)
  
  tmp_genes_to_label <- data_tmp %>% 
    filter(signif_DE == "UP") %>%
    arrange(padj) %>%
    head(head_number) %>%
    pull(mgi_symbol)
  genes_to_label <- c(genes_to_label, tmp_genes_to_label)
  
  tmp_genes_to_label <- data_tmp %>% 
    filter(padj <= padj_cutoff & abs(log2FoldChange) >= log2FC_cutoff) %>%
    arrange(log2FoldChange) %>%
    head(head_number) %>%
    pull(mgi_symbol)
  genes_to_label <- c(genes_to_label, tmp_genes_to_label)
  
  tmp_genes_to_label <- data_tmp %>% 
    filter(padj <= padj_cutoff & abs(log2FoldChange) >= log2FC_cutoff) %>%
    arrange(log2FoldChange) %>%
    tail(head_number) %>%
    pull(mgi_symbol)
  genes_to_label <- c(genes_to_label, tmp_genes_to_label)
  
  genes_to_label <- unique(genes_to_label)
  
  logical_filter <- genes_to_label %in% (data_tmp %>% filter(signif_DE %in% c("DOWN", "UP")) %>% pull(mgi_symbol))
  genes_to_label <- genes_to_label[logical_filter]
  data_tmp$mgi_symbol <- data_tmp$mgi_symbol
  
  titel <- stringr::str_extract(path_file, pattern = "//.*")
  
  if(titel == "deg_sample_name_used_Tumor_WT_NK_TP1_D_vs_Tumor_only_TP1_D.xlsx"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      scale_y_break(c(20, 21), scales = 0.33)
  }else if (titel == "//deg_sample_name_used_Tumor_Ifng_KO_NK_TP1_A_vs_Tumor_only_TP1_A.xlsx"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      ylim(-0.1, 37) +
      scale_x_continuous(breaks = c(-4, 0, 4, 8, 12), limits = c(-6, 12.5))
  }else if (titel == "//deg_sample_name_used_Tumor_Ifng_KO_NK_TP1_D_vs_Tumor_only_TP1_D.xlsx"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      ylim(-0.1, 25) +
      scale_x_continuous(breaks = c(-2.5, 0, 2.5, 5), limits = c(-2.7, 6.4))
  }else if (titel == "//deg_sample_name_used_Tumor_Ifng_KO_NK_TP2_A_vs_Tumor_only_TP2_A.xlsx"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      ylim(-0.1, 33) +
      scale_x_continuous(breaks = c(-6, -3, 0, 3, 6), limits = c(-6.4, 6.1))
  }else if (titel == "//deg_sample_name_used_Tumor_Ifng_KO_NK_TP2_D_vs_Tumor_only_TP2_D.xlsx"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      ylim(-0.1, 16) +
      scale_x_continuous(breaks = c(-5, 0, 5, 10, 15), limits = c(-8, 16))
  }else{
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8))
  }
  
  print(p)
  
  ggsave(filename = paste0("~/workspace/results/review_additions", titel, ".eps"),
         plot = p,
         device="eps"
  )
}
dev.off()


# investigating ENSMUSG00000091228

tmp <- counts(dds_subset_filt, normalized = TRUE)["ENSMUSG00000091228",]
names(tmp) <- dds_subset_filt$full_sample_name
names(tmp) <- dds_subset_filt$condition

tmp




###############
abs_filt_samples=3
padj_cutoff = 0.05
log2FC_cutoff = 0.58 #(FC=1.5); log2FC=1.0 # (FC=2)
var_expl_needed <- 0.6   

if(exists("deg_design")) {rm(deg_design)}

dds_subset <- dds_filt[ , dds_filt$original_experiment %in% "EXP 9.11"]  
dds_subset <- dds_subset[ , !dds_subset$condition %in% "Tumor+Prf1_KO_NK"] 
dds_subset$condition <- str_replace_all(dds_subset$condition, pattern = "\\+", replacement = "_")
dds_subset <- dds_subset[ , !dds_subset$timepoint_cell_harvesting %in% "TP0"] 

#deg_design <- as.formula("~ cell_line_label * timepoint_cell_harvesting + condition")

# dds_subset$sample_name_used <- factor(dds_subset$sample_name_used)
# groups_used <- dds_subset$sample_name_used
# design <- model.matrix(~ 0 + groups_used)
# colnames(design) <- levels(groups_used)

# dds_subset$condition <- factor(dds_subset$condition, levels = c("Tumor_only", "Tumor+Ifng_KO_NK", "Tumor+WT_NK"))
# table(dds_subset$condition)
# dds_subset$cell_line_label <- factor(dds_subset$cell_line_label , levels = c("A", "D"))
# table(dds_subset$cell_line_label)
# dds_subset$timepoint_cell_harvesting <- factor(dds_subset$timepoint_cell_harvesting , levels = c("TP1", "TP2"))
# table(dds_subset$timepoint_cell_harvesting)


dds_subset$sample_name_used <- str_replace_all(dds_subset$sample_name_used, pattern = "\\+", replacement = "_")
dds_subset$sample_name_used <- factor(dds_subset$sample_name_used)
table(dds_subset$sample_name_used)

deg_design <- as.formula("~ sample_name_used")
design(dds_tmp_filt) <- deg_design

compare_these <- list(c("Tumor_Ifng_KO_NK_TP1_A", "Tumor_only_TP1_A"),
                      c("Tumor_Ifng_KO_NK_TP1_D", "Tumor_only_TP1_D"),
                      c("Tumor_Ifng_KO_NK_TP2_A", "Tumor_only_TP2_A"),
                      c("Tumor_Ifng_KO_NK_TP2_D", "Tumor_only_TP2_D"),
                      c("Tumor_WT_NK_TP1_A", "Tumor_only_TP1_A"),
                      c("Tumor_WT_NK_TP1_D", "Tumor_only_TP1_D"),
                      c("Tumor_WT_NK_TP2_A", "Tumor_only_TP2_A"),
                      c("Tumor_WT_NK_TP2_D", "Tumor_only_TP2_D"))

for(comparison_name in compare_these){
  print(comparison_name)
  
  dds_tmp <- dds_subset[ , dds_subset$sample_name_used %in% comparison_name]  
  dds_tmp$sample_name_used <- factor(dds_tmp$sample_name_used, levels = comparison_name[c(2,1)])
  table(dds_tmp$sample_name_used)
  
  dds_tmp_filt <- filterDatasets(dds_tmp, 
                                    abs_filt = TRUE, 
                                    abs_filt_samples = abs_filt_samples)
  design(dds_tmp_filt) <- deg_design
  
  dds_tmp_filt <- DESeq2::estimateSizeFactors(dds_tmp_filt)
  dds_tmp_filt <- DESeq2::DESeq(dds_tmp_filt) # do not replace outliers based on replicates
  
  log2_norm_dds_tmp_filt <- DESeq2::normTransform(dds_tmp_filt)
  vsd_dds_tmp_filt <- DESeq2::vst(dds_tmp_filt, blind = FALSE) 
  
  coeff_name <- resultsNames(dds_tmp_filt)[2]
  dds_tmp_filt_results <- generateResults_upd(dds_object = dds_tmp_filt, 
                                              coeff_name = coeff_name,
                                              cond_numerator = comparison_name[1], 
                                              cond_denominator = comparison_name[2],
                                              cond_variable = "sample_name_used")
  
  openxlsx::write.xlsx(deg_condition_KO_Ly6a_vs_WT_Ly6a_results, 
                       file = file.path("~/workspace/results/review_additions", paste0("deg_", coeff_name, ".xlsx")))
  
}








# =================================================================================================
# declare options 
options(ggrepel.max.overlaps = Inf)
padj_cutoff <- 0.05
log2FC_cutoff <- 0.58

# actual processign
path_data_folder <- "~/workspace/datasets/Plots to edit/" 
path_data_files <- list.files(path = path_data_folder,
                              pattern = "*differential.tsv",
                              full.names = TRUE)
head_number <- 10


pdf(file = "~/workspace/results/review_additions/volcanos_10.pdf")
for(path_file in path_data_files){
  #path_file = "/home/rstudio/workspace/datasets/Plots to edit//balanced_woo_contrast_treatment_IfngKONK_vs_untreated_TP1 A_differential.tsv"
  #path_file = "/home/rstudio/workspace/datasets/Plots to edit//balanced_woo_contrast_treatment_IfngKONK_vs_untreated_TP1 D_differential.tsv"
  #path_file = "/home/rstudio/workspace/datasets/Plots to edit//balanced_woo_contrast_treatment_IfngKONK_vs_untreated_TP2 A_differential.tsv"
  #path_file = "/home/rstudio/workspace/datasets/Plots to edit//balanced_woo_contrast_treatment_IfngKONK_vs_untreated_TP2 D_differential.tsv"
  #path_file = "/home/rstudio/workspace/datasets/Plots to edit//balanced_woo_contrast_treatment_IfngKONK_vs_untreated_TP1 D_differential.tsv"
  
  print(path_file)
  data_tmp <- read.table(file = path_file, header = TRUE, sep = "\t")
  data_tmp <- data_tmp %>% 
    dplyr::filter(!is.na(padj))
  data_tmp$signif_DE <- "NO"
  data_tmp$signif_DE[data_tmp$log2FoldChange > 
                                         log2FC_cutoff & data_tmp$padj < padj_cutoff] <- "UP"
  data_tmp$signif_DE[data_tmp$log2FoldChange < 
                                         -log2FC_cutoff & data_tmp$padj < padj_cutoff] <- "DOWN"
  
  genes_to_label <- c("Ly6a", "Plaat3", "Lgals3bp", "Stat1", "Ccl5")
  
  tmp_genes_to_label <- data_tmp %>% 
    filter(signif_DE == "DOWN") %>%
    arrange(padj) %>%
    head(head_number) %>%
    pull(gene_name)
  genes_to_label <- c(genes_to_label, tmp_genes_to_label)

  tmp_genes_to_label <- data_tmp %>% 
    filter(signif_DE == "UP") %>%
    arrange(padj) %>%
    head(head_number) %>%
    pull(gene_name)
  genes_to_label <- c(genes_to_label, tmp_genes_to_label)
  
  tmp_genes_to_label <- data_tmp %>% 
    filter(padj <= padj_cutoff & abs(log2FoldChange) >= log2FC_cutoff) %>%
    arrange(log2FoldChange) %>%
    head(head_number) %>%
    pull(gene_name)
  genes_to_label <- c(genes_to_label, tmp_genes_to_label)
  
  tmp_genes_to_label <- data_tmp %>% 
    filter(padj <= padj_cutoff & abs(log2FoldChange) >= log2FC_cutoff) %>%
    arrange(log2FoldChange) %>%
    tail(head_number) %>%
    pull(gene_name)
  genes_to_label <- c(genes_to_label, tmp_genes_to_label)
  
  genes_to_label <- unique(genes_to_label)
  
  logical_filter <- genes_to_label %in% (data_tmp %>% filter(signif_DE %in% c("DOWN", "UP")) %>% pull(gene_name))
  genes_to_label <- genes_to_label[logical_filter]
  data_tmp$mgi_symbol <- data_tmp$gene_name
  
  titel <- stringr::str_extract(path_file, pattern = "//.*")
  
  if(titel == "//balanced_woo_WTNK_versus_IfngKONK_for_D_TP1_differential.tsv"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      scale_y_break(c(20, 21), scales = 0.33)
  }else if(titel == "//balanced_woo_contrast_treatment_WTNK_vs_untreated_TP1 D_differential.tsv"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      scale_y_break(c(20, 21), scales = 0.33)
  }else if (titel == "//balanced_woo_contrast_treatment_IfngKONK_vs_untreated_TP1 A_differential.tsv"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      ylim(-0.1, 37) +
      scale_x_continuous(breaks = c(-4, 0, 4, 8, 12), limits = c(-6, 12.5))
  }else if (titel == "//balanced_woo_contrast_treatment_IfngKONK_vs_untreated_TP1 D_differential.tsv"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      ylim(-0.1, 20) +
      scale_x_continuous(breaks = c(-2.5, 0, 2.5, 5), limits = c(-2.7, 6.4))
  }else if (titel == "//balanced_woo_contrast_treatment_IfngKONK_vs_untreated_TP2 A_differential.tsv"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      ylim(-0.1, 33) +
      scale_x_continuous(breaks = c(-6, -3, 0, 3, 6), limits = c(-6.4, 6.1))
  }else if (titel == "//balanced_woo_contrast_treatment_IfngKONK_vs_untreated_TP2 D_differential.tsv"){
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8)) +
      ylim(-0.1, 16) +
      scale_x_continuous(breaks = c(-5, 0, 5, 10, 15), limits = c(-8, 16))
  }else{
    p <- plotVolcano(data_tmp, genes_of_interest = genes_to_label, plot_title = titel) + 
      theme(plot.title = element_text(size = 8))
  }
  
  print(p)
  
  ggsave(filename = paste0("~/workspace/results/review_additions", titel, ".eps"),
         plot = p,
         device="eps"
  )
}
dev.off()


# Enrichments
pdf(file = "~/workspace/results/review_additions/Ly6AE_enrichments.pdf")
path_hyp_export <- "~/workspace/results/review_additions/Ly6AE_enrichments_tables/"

# HypeR GSEA plots

#gs_hallmark <- hypeR::msigdb_gsets(species = "Mus musculus", category = c("H"), clean=TRUE) 
#gs_C2_kegg <- hypeR::msigdb_gsets(species = "Mus musculus", category = c("C2"), subcategory = "CP:KEGG", clean=TRUE) 
gs_C5_GOBP <- hypeR::msigdb_gsets(species = "Mus musculus", category = c("C5"), subcategory = "GO:BP", clean=TRUE) 
gs_C5_GOCC <- hypeR::msigdb_gsets(species = "Mus musculus", category = c("C5"), subcategory = "GO:CC", clean=TRUE) 
#gs_C5_GOMF <- hypeR::msigdb_gsets(species = "Mus musculus", category = c("C5"), subcategory = "GO:MF", clean=TRUE) 

path_Ly6a_KO_vs_WT <- path_data_files[10]
Ly6a_KO_vs_WT <- read.table(file = path_Ly6a_KO_vs_WT,
                            header = TRUE,
                            sep = "\t")
Ly6a_KO_vs_WT <- Ly6a_KO_vs_WT %>% filter(padj <= padj_cutoff,
                                          abs(log2FoldChange) >= 0.58)
Ly6a_KO_vs_WT_up <- Ly6a_KO_vs_WT %>% filter(log2FoldChange > 0)
Ly6a_KO_vs_WT_down <- Ly6a_KO_vs_WT %>% filter(log2FoldChange < 0)

enrichment_C5_GOBP_up <- hypeR::hypeR(signature = Ly6a_KO_vs_WT_up$gene_name, 
                                                         genesets = gs_C5_GOBP, 
                                                         test = "hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOBP_up,
                    file_path = paste0(path_hyp_export, "GOBP_Ly6a_KO_vs_WT_up.xlsx"))
enrichment_C5_GOBP_plot_up <- hyp_dots(enrichment_C5_GOBP_up, 
                                           merge = TRUE, 
                                           fdr = 0.05, 
                                           top = 20, 
                                           abrv = 70, 
                                           val = "fdr", 
                                           title = "GOBP Ly6a_KO_vs_WT_up") +
  theme_bw()
enrichment_C5_GOBP_plot_up
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOBP_Ly6a_KO_vs_WT_up.eps"),
       plot = enrichment_C5_GOBP_plot_up,
       device="eps"
)

enrichment_C5_GOBP_down <- hypeR::hypeR(signature = Ly6a_KO_vs_WT_down$gene_name, 
                                                         genesets = gs_C5_GOBP, 
                                                         test = "hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOBP_down,
                    file_path = paste0(path_hyp_export, "GOBP_Ly6a_KO_vs_WT_down.xlsx"))
enrichment_C5_GOBP_plot_down <- hyp_dots(enrichment_C5_GOBP_down, 
                                           merge = TRUE, 
                                           fdr = 0.05, 
                                           top = 20, 
                                           abrv = 70, 
                                           val = "fdr", 
                                           title = "GOBP Ly6a_KO_vs_WT_down") +
  theme_bw()
enrichment_C5_GOBP_plot_down
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOBP_Ly6a_KO_vs_WT_down.eps"),
       plot = enrichment_C5_GOBP_plot_down,
       device="eps"
)

enrichment_C5_GOCC_up <- hypeR::hypeR(signature = Ly6a_KO_vs_WT_up$gene_name, 
                                        genesets = gs_C5_GOCC, 
                                        test = "hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOCC_up,
                    file_path = paste0(path_hyp_export, "GOCC_Ly6a_KO_vs_WT_up.xlsx"))
enrichment_C5_GOCC_plot_up <- hyp_dots(enrichment_C5_GOCC_up, 
                                         merge = TRUE, 
                                         fdr = 0.05, 
                                         top = 20, 
                                         abrv = 70, 
                                         val = "fdr", 
                                         title = "GOCC Ly6a_KO_vs_WT_up") +
  theme_bw()
enrichment_C5_GOCC_plot_up
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOCC_Ly6a_KO_vs_WT_up.eps"),
       plot = enrichment_C5_GOCC_plot_up,
       device="eps"
)

enrichment_C5_GOCC_down <- hypeR::hypeR(signature = Ly6a_KO_vs_WT_down$gene_name, 
                                        genesets = gs_C5_GOCC, 
                                        test = "hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOCC_down,
                    file_path = paste0(path_hyp_export, "GOCC_Ly6a_KO_vs_WT_down.xlsx"))
enrichment_C5_GOCC_plot_down <- hyp_dots(enrichment_C5_GOCC_down, 
                                         merge = TRUE, 
                                         fdr = 0.05, 
                                         top = 20, 
                                         abrv = 70, 
                                         val = "fdr", 
                                         title = "GOCC Ly6a_KO_vs_WT_down") +
  theme_bw()
enrichment_C5_GOCC_plot_down
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOCC_Ly6a_KO_vs_WT_down.eps"),
       plot = enrichment_C5_GOCC_plot_down,
       device="eps"
)


#gs_hallmark_h <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("H"), clean=TRUE) 
#gs_C2_kegg_h <-  hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:KEGG", clean=TRUE) 
gs_C5_GOBP_h <-  hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:BP", clean=TRUE) 
gs_C5_GOCC_h <-  hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:CC", clean=TRUE) 
#gs_C5_GOMF_h <-  hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:MF", clean=TRUE) 


path_Ly6e_KO_vs_WT <- path_data_files[9]
Ly6e_KO_vs_WT <- read.table(file = path_Ly6e_KO_vs_WT,
                            header = TRUE,
                            sep = "\t")
Ly6e_KO_vs_WT <- Ly6e_KO_vs_WT %>% filter(padj <= padj_cutoff, abs(log2FoldChange) >= 0.58)

Ly6e_KO_vs_WT_up <- Ly6e_KO_vs_WT %>% filter(log2FoldChange > 0)
Ly6e_KO_vs_WT_down <- Ly6e_KO_vs_WT %>% filter(log2FoldChange < 0)

enrichment_C5_GOBP_up <- hypeR::hypeR(signature = Ly6e_KO_vs_WT_up$gene_name, 
                                   genesets = gs_C5_GOBP_h, 
                                   test = "hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOBP_up,
                    file_path = paste0(path_hyp_export, "GOBP_Ly6e_KO_vs_WT_up.xlsx"))
enrichment_C5_GOBP_plot_up <- hyp_dots(enrichment_C5_GOBP_up, 
                                           merge = TRUE, 
                                           fdr = 0.05, 
                                           top = 20, 
                                           abrv = 70, 
                                           val = "fdr", 
                                           title = "GOBP Ly6e_KO_vs_WT_up") +
  theme_bw()
enrichment_C5_GOBP_plot_up
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOBP_Ly6e_KO_vs_WT_up.eps"),
       plot = enrichment_C5_GOBP_plot_up,
       device="eps"
)

enrichment_C5_GOBP_down <- hypeR::hypeR(signature = Ly6e_KO_vs_WT_down$gene_name, 
                                      genesets = gs_C5_GOBP_h, 
                                      test = "hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOBP_down,
                    file_path = paste0(path_hyp_export, "GOBP_Ly6e_KO_vs_WT_down.xlsx"))
enrichment_C5_GOBP_plot_down <- hyp_dots(enrichment_C5_GOBP_down, 
                                       merge = TRUE, 
                                       fdr = 0.05, 
                                       top = 20, 
                                       abrv = 70, 
                                       val = "fdr", 
                                       title = "GOBP Ly6e_KO_vs_WT_down") +
  theme_bw()
enrichment_C5_GOBP_plot_down
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOBP_Ly6e_KO_vs_WT_down.eps"),
       plot = enrichment_C5_GOBP_plot_down,
       device="eps"
)


enrichment_C5_GOCC_up <- hypeR::hypeR(signature = Ly6e_KO_vs_WT_up$gene_name, 
                                      genesets = gs_C5_GOCC_h, 
                                      test = "hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOCC_up,
                    file_path = paste0(path_hyp_export, "GOCC_Ly6e_KO_vs_WT_up.xlsx"))
enrichment_C5_GOCC_plot_up <- hyp_dots(enrichment_C5_GOCC_up, 
                                       merge = TRUE, 
                                       fdr = 0.05, 
                                       top = 20, 
                                       abrv = 70, 
                                       val = "fdr", 
                                       title = "GOCC Ly6e_KO_vs_WT_up") +
  theme_bw()
enrichment_C5_GOCC_plot_up
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOCC_Ly6e_KO_vs_WT_up.eps"),
       plot = enrichment_C5_GOCC_plot_up,
       device="eps"
)


enrichment_C5_GOCC_down <- hypeR::hypeR(signature = Ly6e_KO_vs_WT_down$gene_name, 
                                        genesets = gs_C5_GOCC_h, 
                                        test = "hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOCC_down,
                    file_path = paste0(path_hyp_export, "GOCC_Ly6e_KO_vs_WT_down.xlsx"))
enrichment_C5_GOCC_plot_down <- hyp_dots(enrichment_C5_GOCC_down, 
                                         merge = TRUE, 
                                         fdr = 0.05, 
                                         top = 20, 
                                         abrv = 70, 
                                         val = "fdr", 
                                         title = "GOCC Ly6e_KO_vs_WT_down") +
  theme_bw()
enrichment_C5_GOCC_plot_down
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOCC_Ly6e_KO_vs_WT_down.eps"),
       plot = enrichment_C5_GOCC_plot_down,
       device="eps"
)

# check what are the orthologs of mice genes in human and do GOBP
library(Orthology.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

Human_orthologs_ly6a_up <- mapfun(Ly6a_KO_vs_WT_up$gene_name)
Human_orthologs_ly6a_up <- Human_orthologs_ly6a_up$Human_symbol
Human_orthologs_ly6a_up <- Human_orthologs_ly6a_up[!is.na(Human_orthologs_ly6a_up)]

enrichment_C5_GOBP_up <- hypeR::hypeR(signature = Human_orthologs_ly6a_up, 
                                      genesets = gs_C5_GOBP_h, 
                                      test="hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOBP_up,
                    file_path = paste0(path_hyp_export, "GOBP_Ly6a_KO_vs_WT_UP-regulated_HUMAN_ORTHOLOGUES.xlsx"))
enrichment_C5_GOBP_plot_up <- hyp_dots(enrichment_C5_GOBP_up, 
                                       merge=TRUE, 
                                       fdr=0.05, 
                                       top = 20, 
                                       abrv=70, 
                                       val="fdr", 
                                       title="GOBP Ly6a_KO_vs_WT UP-regulated HUMAN ORTHOLOGUES") +
  theme_bw()
enrichment_C5_GOBP_plot_up
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOBP_Ly6a_KO_vs_WT_UP-regulated_HUMAN_ORTHOLOGUES.eps"),
       plot = enrichment_C5_GOBP_plot_up,
       device="eps"
)

Human_orthologs_ly6a_down <- mapfun(Ly6a_KO_vs_WT_down$gene_name)
Human_orthologs_ly6a_down <- Human_orthologs_ly6a_down$Human_symbol
Human_orthologs_ly6a_down <- Human_orthologs_ly6a_down[!is.na(Human_orthologs_ly6a_down)]

enrichment_C5_GOBP_down <- hypeR::hypeR(signature = Human_orthologs_ly6a_down, 
                                        genesets = gs_C5_GOBP_h, 
                                        test="hypergeometric")
hypeR::hyp_to_excel(hyp_obj = enrichment_C5_GOBP_down,
                    file_path = paste0(path_hyp_export, "GOBP_Ly6a_KO_vs_WT_DOWN-regulated_HUMAN_ORTHOLOGUES.xlsx"))
enrichment_C5_GOBP_plot_down <- hyp_dots(enrichment_C5_GOBP_down, 
                                         merge=TRUE, 
                                         fdr=0.05, 
                                         top = 20, 
                                         abrv=70, 
                                         val="fdr", 
                                         title="GOBP Ly6a_KO_vs_WT down-regulated HUMAN ORTHOLOGUES") +
  theme_bw()
enrichment_C5_GOBP_plot_down
ggsave(filename = paste0("~/workspace/results/review_additions/Ly6AE_plots/enrichment_GOBP_Ly6a_KO_vs_WT_DOWN-regulated_HUMAN_ORTHOLOGUES.eps"),
       plot = enrichment_C5_GOBP_plot_down,
       device="eps"
)

dev.off()






