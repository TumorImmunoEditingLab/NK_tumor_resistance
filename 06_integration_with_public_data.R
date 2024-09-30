import::from(.from = readr, read_csv, cols)
import::from(magrittr, "%>%")
import::from(dplyr, mutate, select, filter, rename, arrange, desc, group_by, summarise, pull, ungroup)  # dplyr_mutate = mutate
import::from(purrr, map)
import::from(future, plan, multisession, sequential)
import::from(furrr, furrr_options, future_map2)
import::from(ggplot2, .all=TRUE) # importing all as there is too many
import::from(grid, gpar) # needed in complexheatmap
import::from(kableExtra, kable_styling, kbl)

import::from(.from = GenomicFeatures, makeTxDbFromGFF)
import::from(.from = AnnotationDbi, annot_db_keys = keys, annot_db_select = select)
import::from(.from = DESeq2, .all=TRUE)
import::from(.from = tximport, tximport)

import::from(.from = here::here("utils/filterDatasets.R"), "filterDatasets", .character_only=TRUE) # used for filtering
import::from(.from = here::here("utils/generateEnsemblAnnotation.R"), "generateEnsemblAnnotation", .character_only=TRUE) # used for filtering
import::from(.from = here::here("utils/generateResults.R"), "meanExprsPerGroup", "generateResults_upd", .character_only=TRUE) # used for filtering


# load DE tables from scRNA-seq work.
# Downloaded from here:
# https://immunogenomics.shinyapps.io/nkheme/
untar("~/workspace/datasets/public_data.tar.gz", exdir = "~/workspace/datasets/")

E697 <- read.delim("~/workspace/datasets/public_data/697_Expanded NK-treated.txt")
kasumi2 <- read.delim("~/workspace/datasets/public_data/KASUMI2_Expanded NK-treated.txt")
nalm6 <- read.delim("~/workspace/datasets/public_data/NALM6_Expanded NK-treated.txt")
rchacv <- read.delim("~/workspace/datasets/public_data/RCHACV_Expanded NK-treated.txt")

# Load human to Mouse gene translation dictionarry
H2m <- read.delim("~/workspace/datasets/public_data/H2M.txt", header = FALSE)  
H2m <- H2m[, 1:4]
colnames(H2m) <- c("H_symbol", "H_entrezId", "M_symbol", "M_entrezId")

# Load our data
load(file.path("~/workspace/RDSs/nk_tum_immunoedit_complete_dds_objects.RData"))
load(file.path("~/workspace/RDSs/nk_tum_immunoedit_complete_log2_vsd_filt_objects.RData"))

# add filtration parameters
abs_filt_samples=3
padj_cutoff = 0.05
log2FC_cutoff = 0.58 #(FC=1.5); log2FC=1.0 # (FC=2)
var_expl_needed <- 0.6         # at least 60% variance explained needed

# compare Tumor_only between timepoints!
# condition_tp_subset <- c("Tumor_only_timepoint_0", 
#                          "Tumor_only_timepoint_1", 
#                          "Tumor_only_timepoint_2", 
#                          "Tumor_plus_NK_timepoint_1", 
#                          "Tumor_plus_NK_timepoint_2",
#                          "Tumor_plus_WT_NK_timepoint_2")

condition_tp_subset <- c("Tumor_only_timepoint_1",
                         "Tumor_plus_NK_timepoint_2",
                         "Tumor_plus_NK_timepoint_1",
                         "Tumor_plus_WT_NK_timepoint_2")


if(exists("deg_design")) {rm(deg_design)}
deg_design <- as.formula("~ experiment + cell_line_label + condition_tp")
deg_name <- "TpNK_tp2_vs_Tonly_tp1"
deg_dir <- file.path(output_dir, paste0("deg_", deg_name, "/"))
dir.create(deg_dir)

experiment_subset <- c("EXP9.4", "EXP9.5", "EXP9.6", "EXP9.7")

# just in case remove any existing subsets
if(exists("dds_subset")) {rm(dds_subset)}
if(exists("dds_subset_filt")) {rm(dds_subset_filt)}
if(exists("vsd_subset")) {rm(vsd_subset)}

dds_subset <- dds_filt[ , dds_filt$condition_tp %in% condition_tp_subset]  
dds_subset <- dds_subset[ , dds_subset$experiment %in% experiment_subset]

# select one of these:


dds_subset <- dds_subset ; out_folder <- Homologies_overlap_ABCD
#dds_subset <- dds_subset[ , dds_subset$cell_line_label %in% c("A", "B")] ; out_folder <- Homologies_overlap_AB
#dds_subset <- dds_subset[ , dds_subset$cell_line_label %in% c("C", "D")] ; out_folder <- Homologies_overlap_CD

# fixing Tumor_plus_WT_NK_timepoint_2 -> Tumor_plus_NK_timepoint_2
dds_subset$condition_tp <- droplevels(dds_subset$condition_tp)
table(dds_subset$condition_tp)

dds_subset$experiment <- droplevels(dds_subset$experiment)
table(dds_subset$experiment)

dds_subset$cell_line_label <- droplevels(dds_subset$cell_line_label)
table(dds_subset$cell_line_label)

# filtering lowly expressed genes
dds_subset_filt <- filterDatasets(
  dds_subset,
  abs_filt = TRUE,
  abs_filt_samples = abs_filt_samples
) # at least in N samples, which is a smallest group size

cat("... dds\n")
design(dds_subset_filt) <- deg_design
dds_subset_filt <- DESeq2::estimateSizeFactors(dds_subset_filt)
dds_subset_filt <- DESeq2::DESeq(dds_subset_filt) # do not replace outliers based on replicates

log2_norm_subset_filt <- DESeq2::normTransform(dds_subset_filt)
vsd_subset_filt <- DESeq2::vst(dds_subset_filt, blind = FALSE) # using blind=FALSE utilize design info; blind = TRUE for QC
#rld_filt <- DESeq2::rlog(dds_subset_filt, blind = FALSE) # more robust to differences in sequencing depth!
#rld_filt - too big for >100 samples; skipping! 


deg_TpNK_tp2_vs_Tonly_tp1_results <- generateResults_upd(
  dds_object = dds_subset_filt,
  coeff_name = "condition_tp_Tumor_plus_NK_timepoint_2_vs_Tumor_only_timepoint_1",
  cond_numerator = "Tumor_plus_NK_timepoint_2", 
  cond_denominator = "Tumor_only_timepoint_1",
  cond_variable = "condition_tp"
)


# Now we have our DE table from AB line and we need to translate Mouse genes to Human ones
match(E697$gene, table = H2m$H_symbol)
H2m$M_entrezId <- stringr::str_replace_all(H2m$M_entrezId, pattern = "MGI:", replacement = "")

deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$Hum_homlg_gene <- 
  H2m$H_symbol[match(deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$mgi_symbol, H2m$M_symbol)]

# mart_import <- read.delim("~/workspace/datasets/public_data/mart_export.tsv")
# mart_import <- mart_import %>% filter(Human.gene.name != "")

# deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$Hum_homlg_gene_mart <- 
#   mart_import$Human.gene.name[match(deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$ensembl_id, mart_import$Gene.stable.ID)]


# mart_import_reduced <- mart_import[, c(5,2)]
# H2m_import_reduced <- H2m[, c(3,1)]
# 
# names(mart_import_reduced) <- c("Mouse", "Human")
# names(H2m_import_reduced) <- c("Mouse", "Human")

# combined_h2m <- rbind(mart_import_reduced, H2m_import_reduced)
# combined_h2m <- combined_h2m[combined_h2m$Human != "", ]
# combined_h2m <- combined_h2m[!is.na(combined_h2m$Human), ]
# dim(unique(combined_h2m))
# 
# deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$Hum_homlg_combined <- 
#   combined_h2m$Human[match(deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$mgi_symbol, unique(combined_h2m)$Mouse)]
# 
#dim(combined_h2m)

intersect_E697 <- E697 %>% filter(., p_adj <= 0.05) %>% pull(gene) %>% intersect(., unique(deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$Hum_homlg_gene))
df_tmp <- data.frame("Mouse_" = deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif %>%
                       filter(., Hum_homlg_gene %in% intersect_E697) %>% 
                       distinct(., Hum_homlg_gene, .keep_all = T) %>% 
                       arrange(., Hum_homlg_gene) %>% 
                       select(Hum_homlg_gene, log2FoldChange, padj, mgi_symbol),
                     "E697_" = E697 %>% filter(., p_adj <= 0.05) %>% 
                       filter(., gene %in% intersect_E697) %>% 
                       distinct(., gene, .keep_all = T) %>% 
                       arrange(., gene) %>%
                       select(avg_log2FC, p_adj, gene)
)
write.table(df_tmp, 
            file = paste("~/workspace/results/", out_folder, "/E697.tsv"),
            sep = "\t",append = F, row.names = F, col.names = T, quote = F)
a1 <- length(intersect_E697)
a2 <- length(unique(deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$Hum_homlg_gene))-1 - a1 #we exclude NA/empty from counting
b1 <- length(E697 %>%  filter(., p_adj <= 0.05) %>% pull(gene)) - a1    
b2 <- 20000 - a1 - a2 - b1
c_table <- data.frame("DE" = c(a1, a2),
                      "non_DE"   = c(b1, b2),
                      row.names = c("in", "out"))
fisher.test(c_table, alternative = "greater")



intersect_kasumi2 <- kasumi2 %>%  filter(., p_adj <= 0.05) %>% pull(gene) %>% intersect(., unique(deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$Hum_homlg_gene))
df_tmp <- data.frame("Mouse_" = deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif %>%
                       filter(., Hum_homlg_gene %in% intersect_kasumi2) %>% 
                       distinct(., Hum_homlg_gene, .keep_all = T) %>% 
                       arrange(., Hum_homlg_gene) %>% 
                       select(Hum_homlg_gene, log2FoldChange, padj, mgi_symbol),
                     "kasumi2_" = kasumi2 %>% filter(., p_adj <= 0.05) %>% 
                       filter(., gene %in% intersect_kasumi2) %>% 
                       distinct(., gene, .keep_all = T) %>% 
                       arrange(., gene) %>%
                       select(avg_log2FC, p_adj, gene)
)
write.table(df_tmp, 
            file = paste("~/workspace/results/", out_folder, "/E697.tsv"),
            sep = "\t",append = F, row.names = F, col.names = T, quote = F)

a1 <- length(intersect_kasumi2)
a2 <- length(unique(deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$Hum_homlg_gene))-1 - a1 #we exclude NA/empty from counting
b1 <- length(kasumi2 %>%  filter(., p_adj <= 0.05) %>% pull(gene)) - a1    
b2 <- 20000 - a1 - a2 - b1
c_table <- data.frame("DE" = c(a1, a2),
                      "non_DE" = c(b1, b2),
                      row.names = c("in", "out"))
fisher.test(c_table, alternative = "greater")



intersect_nalm6 <- nalm6 %>%  filter(., p_adj <= 0.05) %>% pull(gene) %>% intersect(.,deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$Hum_homlg_gene)
df_tmp <- data.frame("Mouse_" = deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif %>%
                       filter(., Hum_homlg_gene %in% intersect_nalm6) %>% 
                       distinct(., Hum_homlg_gene, .keep_all = T) %>% 
                       arrange(., Hum_homlg_gene) %>% 
                       select(Hum_homlg_gene, log2FoldChange, padj, mgi_symbol),
                     "nalm6_" = nalm6 %>% filter(., p_adj <= 0.05) %>% 
                       filter(., gene %in% intersect_nalm6) %>% 
                       distinct(., gene, .keep_all = T) %>% 
                       arrange(., gene) %>%
                       select(avg_log2FC, p_adj, gene)
)
write.table(df_tmp, 
            file = paste("~/workspace/results/", out_folder, "/E697.tsv"),
            sep = "\t",append = F, row.names = F, col.names = T, quote = F)

a1 <- length(intersect_nalm6)
a2 <- length(unique(deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$Hum_homlg_gene))-1 - a1 #we exclude NA/empty from counting
b1 <- length(nalm6 %>%  filter(., p_adj <= 0.05) %>% pull(gene)) - a1    
b2 <- 20000 - a1 - a2 - b1
c_table <- data.frame("DE" = c(a1, a2),
                      "non_DE" = c(b1, b2),
                      row.names = c("in", "out"))
fisher.test(c_table, alternative = "greater")





intersect_rchacv <- rchacv %>% filter(., p_adj <= 0.05) %>% pull(gene) %>% intersect(.,deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$Hum_homlg_gene)
df_tmp <- data.frame("Mouse_" = deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif %>%
                       filter(., Hum_homlg_gene %in% intersect_rchacv) %>% 
                       distinct(., Hum_homlg_gene, .keep_all = T) %>% 
                       arrange(., Hum_homlg_gene) %>% 
                       select(Hum_homlg_gene, log2FoldChange, padj, mgi_symbol),
                     "rchacv_" = rchacv %>% filter(., p_adj <= 0.05) %>% 
                       filter(., gene %in% intersect_rchacv) %>% 
                       distinct(., gene, .keep_all = T) %>% 
                       arrange(., gene) %>%
                       select(avg_log2FC, p_adj, gene)
)
write.table(df_tmp, 
            file = paste("~/workspace/results/", out_folder, "/E697.tsv"),
            sep = "\t",append = F, row.names = F, col.names = T, quote = F)

a1 <- length(intersect_rchacv)
a2 <- length(unique(deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$Hum_homlg_gene))-1 - a1 #we exclude NA/empty from counting
b1 <- length(rchacv %>%  filter(., p_adj <= 0.05) %>% pull(gene)) - a1    
b2 <- 20000 - a1 - a2 - b1
c_table <- data.frame("DE" = c(a1, a2),
                      "non_DE" = c(b1, b2),
                      row.names = c("in", "out"))
fisher.test(c_table, alternative = "greater")



# # using combined table
# E697 %>% filter(., p_adj <= 0.05) %>% pull(gene) %>% intersect(., deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$Hum_homlg_combined)
# kasumi2 %>%  filter(., p_adj <= 0.05) %>% pull(gene) %>% intersect(., deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$Hum_homlg_combined)
# nalm6 %>%  filter(., p_adj <= 0.05) %>% pull(gene) %>% intersect(.,deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$Hum_homlg_combined)
# rchacv %>% filter(., p_adj <= 0.05) %>% pull(gene) %>% intersect(.,deg_TpNK_tp2_vs_Tonly_tp1_results$results_signif$Hum_homlg_combined)






