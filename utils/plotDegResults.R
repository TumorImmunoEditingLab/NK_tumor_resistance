# functions to plot DEG results

plotVolcano <- function(dds_results_obj=NULL, genes_of_interest=NULL, plot_title=NULL){
  require(ggrastr)
  
  results_data_annot_forPlot <- dds_results_obj %>%
    dplyr::filter(!is.na(pvalue))
  results_data_annot_forPlot$signif_DE <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  results_data_annot_forPlot$signif_DE[results_data_annot_forPlot$log2FoldChange > log2FC_cutoff & results_data_annot_forPlot$padj < padj_cutoff] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  results_data_annot_forPlot$signif_DE[results_data_annot_forPlot$log2FoldChange < -log2FC_cutoff & results_data_annot_forPlot$padj < padj_cutoff] <- "DOWN"
  results_data_annot_forPlot$signif_DE <- factor(results_data_annot_forPlot$signif_DE,
                                                 levels = c("NO", "DOWN", "UP"))
  table(results_data_annot_forPlot$signif_DE)
  
  signif_volcanoPlot <- ggplot(data = results_data_annot_forPlot, aes(x = log2FoldChange, y = -log10(padj), col=signif_DE)) +
    ggrastr::geom_point_rast() +
    #gghighlight::gghighlight(signif_DE %in% c("DOWN", "UP")) +
    ggrepel::geom_label_repel(data = . %>% filter(mgi_symbol %in% genes_of_interest), aes(label = mgi_symbol),
                              show.legend = FALSE) +
    #geom_vline(xintercept=c(-log2FC_cutoff, log2FC_cutoff), col="red", linetype="dashed") +
    #geom_hline(yintercept=-log10(padj_cutoff), col="red", linetype="dashed") + # need to adjust to match padj_cutoff
    #scale_color_manual(values=c(DOWN="navy", UP="firebrick3")) +
    scale_color_manual(values=c(DOWN="navy", UP="firebrick3", NO = "grey")) +
    theme_bw(base_size = 14) +
    labs(x = "log2FC") + 
    ggtitle(plot_title)
  #y = "-log10( p-value )",color = "signif. DE") 
}


plotVolcano_repel <- function(dds_results_obj=NULL, genes_of_interest=NULL, genes_of_interest2 = NULL, plot_title=NULL, max_overlaps = 20){
  results_data_annot_forPlot <- dds_results_obj %>%
    dplyr::filter(!is.na(pvalue))
  results_data_annot_forPlot$signif_DE <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  results_data_annot_forPlot$signif_DE[results_data_annot_forPlot$log2FoldChange > log2FC_cutoff & results_data_annot_forPlot$padj < padj_cutoff] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  results_data_annot_forPlot$signif_DE[results_data_annot_forPlot$log2FoldChange < -log2FC_cutoff & results_data_annot_forPlot$padj < padj_cutoff] <- "DOWN"
  
  results_data_annot_forPlot$signif_DE[(results_data_annot_forPlot$mgi_symbol %in% names(genes_of_interest2)) & (results_data_annot_forPlot$log2FoldChange > log2FC_cutoff & results_data_annot_forPlot$padj < padj_cutoff)] <- "Exp9.8_specific_UP"
  
  results_data_annot_forPlot$signif_DE[(results_data_annot_forPlot$mgi_symbol %in% names(genes_of_interest2)) & (results_data_annot_forPlot$log2FoldChange < - log2FC_cutoff & results_data_annot_forPlot$padj < padj_cutoff)] <- "Exp9.8_specific_DOWN"
  
  results_data_annot_forPlot$signif_DE <- factor(results_data_annot_forPlot$signif_DE,
                                                 levels = c("NO", "DOWN", "UP", "Exp9.8_specific_UP", "Exp9.8_specific_DOWN"))
  table(results_data_annot_forPlot$signif_DE)
  
  signif_volcanoPlot <- ggplot(data = results_data_annot_forPlot, aes(x = log2FoldChange, y = -log10(pvalue), col=signif_DE)) +
    geom_point() +
    geom_hline(yintercept=-log(padj_cutoff), col="red", linetype="dashed") + # need to adjust to match padj_cutoff
    #gghighlight::gghighlight(signif_DE %in% c("DOWN", "UP")) +
    ggrepel::geom_label_repel(data = . %>% 
                                filter(mgi_symbol %in% c(names(genes_of_interest), names(genes_of_interest2))) %>% 
                                filter(signif_DE %in% c("DOWN", "UP", "Exp9.8_specific_UP", "Exp9.8_specific_DOWN")), 
                                aes(label = mgi_symbol), 
                              show.legend = FALSE,
                              box.padding = 0.5, 
                              segment.color = "black", 
                              min.segment.length = 0,
                              max.overlaps = max_overlaps) +
    scale_color_manual(values = c(DOWN="#553388", UP="#DD3344", NO = "grey", Exp9.8_specific_UP = "#FF8F4A", Exp9.8_specific_DOWN = "#76C9EC")) +
    #geom_vline(xintercept=c(-log2FC_cutoff, log2FC_cutoff), col="red", linetype="dashed") +
    #scale_color_manual(values=c(DOWN="navy", UP="firebrick3")) +
    theme_bw(base_size = 14) +
    labs(x = "log2FC") + 
    ggtitle(plot_title)
  #y = "-log10( p-value )",color = "signif. DE") 
}
