#project setup
project_name <- "DNA_barcode"
base_dir <- "/home/mohamed/data_isilon/Research/KOENIG/zKoenig_Halbritter/"
seq_round <- "EXP 9.11 Revision experiment"
input_dir <- Sys.glob(file.path(base_dir, seq_round, "*","Alignment_1"))

out_dir <- file.path(here::here("new_runs"),paste0("out_",seq_round))
fwr_primer <- "GAACACTCGAGATCAG"
nonvar <- "TGTGGTATGATGT"

#report setup
plots_dir <- glue::glue("{out_dir}/plots/")
docs_dir <- glue::glue("{out_dir}/docs/")
cache_dir <- glue::glue("{out_dir}/cache/")
# dir.create(cache_dir,recursive = TRUE)
# dir.create(plots_dir,recursive = TRUE)
# dir.create(docs_dir,recursive = TRUE)

#load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(stringdist))
suppressPackageStartupMessages(library(VGAM))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ComplexHeatmap))
theme_set(theme_cowplot())


#######
## ---- metadata
meta_data <- file.path(base_dir, seq_round, "EXP9.11_Metadata.csv") %>%
  read_csv() %>%
  janitor::clean_names()%>%
  dplyr::rename(miseq_run = "mi_seq_run_name")%>%
  mutate(sample_name_long = str_replace_all(individual_sample_name_mi_seq_run_experiment_cell_line_a_samplebio_timepoint1_well1_techrep1, " ", ""),
         bio_sample_name = str_replace(sample_name_long, "_R[1|2]", ""),
         cell_line = cell_line_lab_internal_name,
         timepoint = str_replace(timepoint, " ", "_"))

#summary of metadata
metadata_summary_path <- file.path(docs_dir,glue::glue("data_summary_{seq_round}.csv"))
if(!file.exists(metadata_summary_path)){
  meta_data %>%
    dplyr::count(experiment, cell_line, timepoint, sample_bio) %>%
    as.data.frame %>% 
    arrange(timepoint) %>%
    write_csv(metadata_summary_path)
} 

## ---- read fastq files
#sequencing directories
seq_dirs <- input_dir
#sequence_files
seq_files <- list.files(seq_dirs, recursive=TRUE, pattern = "_R1_001.fastq.gz",full.names=TRUE) %>% str_subset("Undetermined", negate = TRUE)
names(seq_files) <- str_split(seq_files, "/",simplify = TRUE)[, ncol(str_split(seq_files, "/",simplify = TRUE))] %>% str_replace("_.+", "") %>% str_replace("-","_")

seq_files <- seq_files[meta_data$sample_name]

#read_files
if (!file.exists(glue("{out_dir}/cache/reads_list_EXP9.11.rds"))) {
  reads_list <- lapply(seq_files, function(file){
    #read the sample 
    readDNAStringSet(file, format = "fastq")
  })
  
  #save read_files cache
  write_rds(reads_list,
            file = glue("{out_dir}/cache/reads_list_EXP9.11.rds"))
  
} else {
  reads_list <- read_rds(glue("{out_dir}/cache/reads_list_EXP9.11.rds"))
}


########
## ----  match_sequence
if (!file.exists(glue("{out_dir}/cache/match_list_EXP9.11.rds"))) {
  match_list <- lapply(reads_list, function(dna_f){
    #location of the forward primer exact match in the reads
    primer_match<- vmatchPattern(fwr_primer,
                                 dna_f,
                                 fixed = TRUE)
    #location of the non-variable region exact match in the reads
    nonvar_match<- vmatchPattern(nonvar,
                                 dna_f,
                                 fixed = TRUE)
    list(primer=primer_match,
         nonvar=nonvar_match)
  })
  
  #save_matching_cache
  write_rds(match_list,
            path = glue("{out_dir}/cache/match_list_EXP9.11.rds"))
} else {
  match_list <- read_rds(glue("{out_dir}/cache/match_list_EXP9.11.rds"))
}

## ---- match_position
if (!file.exists(glue("{docs_dir}/match_all_reads/all_matchings_EXP9.11.csv"))) {
  match_position <- lapply(names(match_list), function(sample_id){
    message(sample_id)
    #convert matching to a dataframe
    read_df <- lapply(match_list[[sample_id]], function(x){
      unlist(x) %>%
        as.data.frame 
    })
    
    #start and end of the read head
    match_pos <- lapply(names(read_df), function(x){
      read_df[[x]] %>% 
        rename_with( ~  glue("{x}_{.}"), .cols = -"names")
    }) %>% 
      purrr::reduce(full_join, by= "names")
    
    match_pos
  })
  names(match_position) <- names(match_list)
  
  ## ---- match_summary
  match_df <- map(names(match_list), function(sample_id){
    message(sample_id)
    
    read_id <- reads_list[[sample_id]]@ranges@NAMES
    
    #matching table 
    match_position[[sample_id]] %>% 
      full_join(data.frame(names = read_id)) %>% 
      mutate(#calculate barcode length
        bc_width = nonvar_start - primer_end -1,
        #barcode sequence
        barcode_seq = substring(reads_list[[sample_id]][names],
                                primer_end+1,
                                primer_end+bc_width),
        #exclusion criteria:
        #1- Either read-head or nonvariable-region are not detected.
        #2- Distance between read-head and nonvariable-region is less than expected barcode length (21).
        #3- Barcode contains unidentified bases (N).
        match_quality = ifelse(is.na(bc_width)| bc_width!= 21,"invalid", ifelse(str_detect(barcode_seq, "N"), "invalid", "valid"))
      ) %>%
      select(names,
             primer_start,
             primer_end,
             nonvar_start,
             nonvar_end,
             barcode_seq,
             bc_width,
             match_quality)
  })
  names(match_df) <- names(match_list)
  
  #save_match_df
  match_df %>%
    bind_rows(.id = "sample_name") %>%
    left_join(meta_data %>% dplyr::select(sample_name, cell_line,sample_name_long)) %>%
    dplyr::select(sample_name_long,  cell_line, everything(),-sample_name) %>%
    vroom::vroom_write(glue("{docs_dir}/match_all_reads/all_matchings_EXP9.11.csv"))
} else {
  match_df <- vroom::vroom(glue("{docs_dir}/match_all_reads/all_matchings_EXP9.11.csv"))
}
## ---- filteration_summary
if (!file.exists(glue("{docs_dir}/match_all_reads/match_summary_EXP9.11.csv"))) {
  match_summary <- match_df %>%
    group_by(sample_name_long, match_quality, bc_width) %>% 
    summarize(n = n_distinct(names)) %>%
    ungroup()
  
  match_summary %>%
    left_join(meta_data %>% dplyr::select(sample_name, cell_line,sample_name_long)) %>%
    dplyr::select(sample_name_long,  cell_line, everything(),-sample_name) %>%
    write_csv(glue("{docs_dir}/match_all_reads/match_summary_EXP9.11.csv"))
}

if (!file.exists(glue("{docs_dir}/match_all_reads/filt_summary_EXP9.11.csv"))) {
  filt_summary <- match_summary %>% 
    group_by(sample_name_long, match_quality) %>% 
    summarize(n = sum(n)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = match_quality,values_from = n) %>% 
    mutate(raw = invalid+valid,
           "valid/raw" = scales::percent(valid/raw, accuracy = 2))
  
  filt_summary %>%
    left_join(meta_data %>% dplyr::select(sample_name, cell_line,sample_name_long)) %>%
    dplyr::select(sample_name_long,  cell_line, everything(),-sample_name) %>%
    write_csv(glue("{docs_dir}/match_all_reads/filt_summary_EXP9.11.csv"))
}

## ---- valid_matchings_count
valid_barcodes_location <- match_df %>%
  filter(match_quality == "valid") %>% 
  dplyr::count(sample_name_long, primer_end, nonvar_start) %>% 
  arrange(sample_name_long, nonvar_start)

if(!file.exists(glue("{docs_dir}/match_all_reads/valid_match_location_summary_EXP9.11.csv"))){
  valid_barcodes_location %>%
    left_join(meta_data %>% dplyr::select(sample_name, cell_line,sample_name_long)) %>%
    dplyr::select(sample_name_long,  cell_line, everything(),-sample_name) %>%
    write_csv(glue("{docs_dir}/match_all_reads/valid_match_location_summary_EXP9.11.csv"))
}

## remove objects
rm(reads_list)          
rm(match_list)          
rm(match_summary)

#plots
pdf(glue("{plots_dir}/qc_valid_reads.pdf"))
lapply(unique(meta_data$experiment), function(cline_exper){
  bind_rows(filt_summary, filt_summary %>% mutate(group = "percent")) %>%
    mutate(group = ifelse(is.na(group), "count", "percent")) %>%
    left_join(meta_data) %>%
    filter(experiment == cline_exper) %>% 
    arrange(cell_line, bio_sample_name) %>%
    mutate(sample_name = factor(sample_name, unique(sample_name))) %>%
    ggplot()+
    geom_col(data = ~subset(.x, group == "count"),aes(sample_name_long, raw, fill = cell_line), alpha = 0.5)+
    geom_col(data = ~subset(.x, group == "count"),aes(sample_name_long, valid, fill = cell_line ))+
    geom_text(data = ~subset(.x, group == "percent"),aes(sample_name_long, valid, label = `valid/raw`),size = 2)+
    scale_color_manual(values = c("13G" = "#d6604d","15O" = "#878787"))+
    coord_flip()+
    labs(x = "Sample name",
         y = "Number of reads",
         title = "Quality control",
         subtitle = "percentage of valid reads")+
    theme(axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 7))
})
dev.off()

############
# filter match
## ----  filter_valid
valid_match <- match_df %>%
  #filter valid reads
  filter(match_quality == "valid")
## ---- count_valid
valid_match_count <- valid_match %>%
  #count barcode
  dplyr::count(sample_name_long,barcode_seq) %>%
  left_join(meta_data %>% dplyr::select(sample_name, cell_line, timepoint, experiment, sample_name_long)) %>%
  dplyr::select(cell_line, timepoint, experiment, sample_name_long, barcode_seq, n) 

rm(valid_match)
rm(match_df)
## ----
sample_lst <- meta_data %>%
  mutate(sample_bio = factor(sample_bio, levels = c("Tumor_only", "Tumor+IFNg KO NK cells", "Tumor+Perforin KO NK cells", "Tumor+NK cells"))) %>%
  arrange(sample_bio,timepoint) %$%
  split(sample_name_long, paste0(cell_line,"_",experiment))  

# aggregate replicates
cline_exper_all_aggr_count <- lapply(names(sample_lst), function(cell_line_expr){
  valid_match_count %>% 
    filter(sample_name_long %in% sample_lst[[cell_line_expr]]) %>%
    dplyr::select(barcode_seq, sample_name_long, n) %>%
    mutate(sample_name_bio = str_replace(sample_name_long, "_R.$", ""),
           sample_name = NULL) %>%
    group_by(barcode_seq, sample_name_bio) %>%
    summarize(barcode_seq = unique(barcode_seq),
              sample_name_bio = unique(sample_name_bio),
              n = sum(n)) %>% #sum barcode counts across the two replicates
    ungroup() %>%
    pivot_wider(names_from ="sample_name_bio", values_from = "n", values_fill = 0)}) %>%
  setNames(., names(sample_lst))

# blue-list barcodes aggregate replicates
cline_exper_blue_aggr_count <- lapply(names(cline_exper_all_aggr_count), function(cline){
  cline_exper_all_aggr_count[[cline]]%>%
    #barcode detected across the 3 tumor only aggregated samples and have a sum of atleast 9 reads
    filter(apply(dplyr::select(., contains("Tumor_only_timepoint0"))>=1,1,sum )==3 & apply(dplyr::select(., contains("Tumor_only_timepoint0")),1,sum )>=9)
}) %>%
  setNames(., names(cline_exper_all_aggr_count))

# filter blue-list barcodes specific to a cell-line, experiment combination
cline_exper_blue_count <- lapply(names(sample_lst), function(cell_line_expr){
  blue_barcodes <- cline_exper_blue_aggr_count[[cell_line_expr]]$barcode_seq
  valid_match_count %>% 
    filter(sample_name_long %in% sample_lst[[cell_line_expr]] & barcode_seq %in% blue_barcodes) %>% #select blue-list barcodes in the cell-line & experiment samples
    dplyr::select(barcode_seq, sample_name_long, n) %>%
    pivot_wider(names_from ="sample_name_long", values_from = "n", values_fill = 0)}) %>%
  setNames(., names(sample_lst))

# barcodes concordance between technical replicates
rep_n <- lapply(names(sample_lst), function(cell_line_expr){
  blue_barcodes <- cline_exper_blue_aggr_count[[cell_line_expr]]$barcode_seq
  
  valid_match_count %>% 
    filter(sample_name_long %in% sample_lst[[cell_line_expr]] ) %>%
    dplyr::select(barcode_seq, sample_name_long, n) %>%
    mutate(reps = str_extract(sample_name_long, "R."),
           sample_name_bio = str_replace(sample_name_long, "_R.$", ""),
           sample_name_long = NULL,
           is_blue = barcode_seq %in% blue_barcodes) %>% 
    pivot_wider(names_from = "reps",
                values_from = "n",
                values_fill = list(n = FALSE)) %>% 
    mutate(detection = apply(.[,c("R1", "R2")], 1, function(x){ ifelse(sum(x>0)==2, "both", names(x)[x>0])}))
}) %>% 
  setNames(names(sample_lst))

pdf(glue("{plots_dir}/technical_reps_scatter.pdf"))
rep_n %>% 
  bind_rows(.id = "group") %>% 
  ggplot()+
  ggrastr::geom_point_rast(aes(R1, R2, color = detection, shape = is_blue),
                           size = 1)+
  facet_wrap(~ group, ncol = 3)+
  scale_shape_manual(values = c(`TRUE`=16,`FALSE`=4))+
  coord_fixed()+  
  guides(colour = guide_legend(override.aes = list(size = 4,
                                                   alpha = 1)),
         shape = guide_legend(override.aes = list(size = 4,
                                                  alpha = 1)))+
  labs(x = "R1 counts",
       y = "R2 counts")

rep_n %>% 
  bind_rows(.id = "group") %>% 
  mutate(R1 = log2(R1+1),
         R2 = log2(R2+1)) %>%
  ggplot()+
  ggrastr::geom_point_rast(aes(R1, R2, color = detection, shape = is_blue),
                           size = 1,
                           alpha = 0.5)+
  facet_wrap(~ group, ncol = 3)+
  scale_shape_manual(values = c(`TRUE`=16,`FALSE`=4))+
  coord_fixed()+
  guides(colour = guide_legend(override.aes = list(size = 4,
                                                   alpha = 1)),
         shape = guide_legend(override.aes = list(size = 4,
                                                  alpha = 1)))+
  labs(x = "log2(R1 counts + 1)",
       y = "log2(R2 counts + 1)")

dev.off()

pdf(glue("{plots_dir}/count_patterns.pdf"), height=8, width = 4)
# Patterns of aggregated counts of blue list barcodes
lapply(names(cline_exper_blue_aggr_count),function(x){
  cols <- structure(1:4, names = c("0","1", "2", ">=3"))
  
  hm_mat <- cline_exper_blue_aggr_count[[x]] %>% 
    dplyr::select(barcode_seq, matches("Tumor_only_timepoint0")) %>% 
    as.data.frame %>% 
    column_to_rownames("barcode_seq") %>% 
    mutate_all(~ifelse(.x>=3,">=3",as.character(.x))) %>%
    group_by_all() %>% 
    dplyr::count() %>% 
    ungroup() %>% 
    arrange_all() %>% 
    arrange(desc(n)) 
  
  ComplexHeatmap::Heatmap(hm_mat[,!colnames(hm_mat)=="n"],
                          cluster_rows=FALSE,
                          cluster_columns=FALSE, 
                          column_title  = x,
                          name = "raw counts",
                          col = cols,
                          right_annotation = rowAnnotation(`log10(n)` = anno_barplot(log10(hm_mat[,"n"]))))    
  
})
dev.off()

# shannon diversity
shannon_diversity_lst <- lapply(cline_exper_blue_aggr_count, function(cline_exper){
  apply(cline_exper[,-1], 2,  function(x) {
    p <- x/sum(x)
    sh_div <- sum(p * log(p), na.rm = TRUE) * -1
    sh_div
  })
})

pdf(glue("{plots_dir}/shannon_diversity.pdf"))
names(shannon_diversity_lst) %>%
  lapply(function(cline_exper){
    shannon_diversity_lst[[cline_exper]] %>%
      enframe %>%
      left_join(meta_data %>% distinct(bio_sample_name, sample_bio), by = c("name" = "bio_sample_name")) %>%
      dplyr::arrange(value) %>%
      mutate(name = fct_inorder(name))%>%
      ggplot()+
      geom_col(aes(name, value, fill = sample_bio))+
      labs(title = cline_exper,
           y = "Shannon diversity",
           x = "",
           fill = "")+
      theme(axis.text.x = element_text(angle = 90),
            legend.position = "top")
  })
dev.off()

## ---- save_tables

valid_match_count %>%
  write_csv(glue("{docs_dir}/match_all_reads/valid_match_count_9.11.csv"))

valid_match_count %>%
  mutate(cline_exp = paste0(cell_line, "-", experiment)) %>% 
  split(.,.$cline_exp) %>% 
  lapply(function(x){
    x %>% 
      dplyr::select(barcode_seq, sample_name_long, n) %>% 
      pivot_wider(names_from ="sample_name_long", values_from = "n", values_fill = 0)
  }) %>% 
  WriteXLS::WriteXLS("{docs_dir}/counts/valid_match_count_matrix_9.11.xlsx")

valid_match_count %>%
  dplyr::select(barcode_seq, sample_name_long, n) %>%
  pivot_wider( names_from = "sample_name_long", values_from = "n") %>%
  mutate_if(is.numeric,coalesce,0) %>%
  write_csv(glue("{docs_dir}/match_all_reads/valid_match_count_matrix_9.11.csv"))

cline_exper_blue_count %>%
  openxlsx::write.xlsx(glue("{docs_dir}/counts/cline_exper_blue_count_9.11.xlsx"))

cline_exper_blue_count %>%
  purrr::reduce(full_join, by = "barcode_seq") %>%
  mutate_if(is.numeric,coalesce,0) %>%
  write_csv(glue("{docs_dir}/counts/blue_list_matrix_9.11.csv"))

cline_exper_blue_aggr_count %>%
  openxlsx::write.xlsx(glue("{docs_dir}/counts/cline_exper_blue_aggr_count_9.11.xlsx"))

cline_exper_blue_aggr_count %>%
  purrr::reduce(full_join, by = "barcode_seq") %>%
  mutate_if(is.numeric,coalesce,0) %>%
  write_csv(glue("{docs_dir}/counts/blue_list_aggr_count_matrix_9.11.csv"))

rep_n %>% 
  bind_rows(.id = "group") %>% 
  write_csv(glue("{docs_dir}/counts/technical_reps_data.9.11.csv"))

rep_n %>% 
  bind_rows(.id = "group") %>% 
  dplyr::count(group,is_blue, detection) %>% 
  pivot_wider(values_from = "n", names_from = "detection")  %>%
  write_csv(glue("{docs_dir}/counts/technical_reps_data_summary.9.11.csv"))

shannon_diversity_lst %>% 
  stack %>%
  dplyr::rename(group = "ind",
                shannon_diversity = "values") %>% 
  rownames_to_column("sample_id") %>% 
  dplyr::select(sample_id, group, shannon_diversity)  %>%
  write_csv(glue("{docs_dir}/counts/blue_aggr_shannon_diversity.9.11.csv")) 



##########
# DEseq2 analysis
#count data
cts_cline_exper <- openxlsx::getSheetNames(glue("{docs_dir}/counts/cline_exper_blue_aggr_count_9.11.xlsx")) %>%
  lapply(function(sheet){
    openxlsx::read.xlsx(glue("{docs_dir}/counts/cline_exper_blue_aggr_count_9.11.xlsx"),
                        sheet = sheet)%>%
      column_to_rownames("barcode_seq")
  })
names(cts_cline_exper) <- openxlsx::getSheetNames(glue("{docs_dir}/counts/cline_exper_blue_aggr_count_9.11.xlsx"))

###deseq metadata
t_order <- c("timepoint_0","timepoint_1", "timepoint_2")
condition_t <- c(paste("Tumor_only",0:2, sep = "_"),
                 "Tumor+NK cells_1",
                 "Tumor+NK cells_2",
                 "Tumor+Perforin KO NK cells_1",
                 "Tumor+Perforin KO NK cells_2",
                 "Tumor+IFNg KO NK cells_1",
                 "Tumor+IFNg KO NK cells_2")
coldata <- meta_data %>%
  mutate(timepoint  = factor(timepoint, levels = t_order),
         sample_bio_t = paste0(sample_bio, "_", str_sub(timepoint,-1)),
         sample_bio_t = factor(sample_bio_t, levels = condition_t)) %>%
  arrange(cell_line, sample_bio_t, timepoint, well) %>%
  dplyr::select(bio_sample_name, sample_bio, timepoint, sample_bio_t, cell_line, well, experiment) %>%
  distinct_all() %>%
  na.omit() %>%
  as.data.frame()
rownames(coldata) <- coldata$bio_sample_name
sample_list <- split(coldata$bio_sample_name, paste0(coldata$cell_line, "_", coldata$experiment))

### comparisons of interest
comp_9.11 <- list(c("Tumor_only_2", "Tumor_only_0"),
                  c("Tumor+NK cells_2", "Tumor_only_0"),
                  c("Tumor+Perforin KO NK cells_2", "Tumor_only_0"),
                  c("Tumor+IFNg KO NK cells_2", "Tumor_only_0"),
                  c("Tumor+NK cells_2", "Tumor_only_2"),
                  c("Tumor+Perforin KO NK cells_2", "Tumor_only_2"),
                  c("Tumor+IFNg KO NK cells_2", "Tumor_only_2"))
comparisons_list <- list(comp_9.11 = comp_9.11)
sample_bio_t_list <- coldata %>% 
  split(paste0(.$cell_line, "_", .$experiment), rownames(.)) %>%
  lapply(function(x){
    split(rownames(x),x$sample_bio_t)
  })

## run deseq for each cell line in each experiment individually
### run deseq
dds_cline_exper <- lapply(names(cts_cline_exper), function(cline){
  cline_exper <- cts_cline_exper[[cline]]
  DESeq2::DESeqDataSetFromMatrix(countData = cline_exper,
                                 colData = coldata[colnames(cline_exper),] %>% mutate(sample_bio_t = droplevels(sample_bio_t)),
                                 design = ~sample_bio_t) %>%
    DESeq2::DESeq()
})

names(dds_cline_exper) <- names(cts_cline_exper)

dds_cline_exper_tumor_ts <- lapply(cts_cline_exper, function(cline_exper){
  DESeq2::DESeqDataSetFromMatrix(countData = cline_exper[, str_subset(colnames(cline_exper), "Tumor_only")],
                                 colData = coldata[str_subset(colnames(cline_exper), "Tumor_only"), ]%>% mutate(sample_bio_t = droplevels(sample_bio_t)),
                                 design = ~sample_bio_t) %>%
    DESeq2::DESeq(test="LRT", reduced=~1)
})               

### results of dds
padj_cutoff <- 0.05
log2fc_cutoff <- 1 
res_lst <- lapply(names(dds_cline_exper), function(cline){
  lapply(comp_9.11, function(comparison){
    DESeq2::results(dds_cline_exper[[cline]], contrast=c("sample_bio_t",comparison))
  })%>%
    setNames(.,sapply(comp_9.11, paste, collapse = "_vs_"))
})%>%
  setNames(.,names(dds_cline_exper))

res_ts_lst <- lapply(names(dds_cline_exper_tumor_ts), function(cline){
  DESeq2::results(dds_cline_exper_tumor_ts[[cline]])
})%>%
  setNames(.,names(dds_cline_exper_tumor_ts))

#add time-series DESeq2 restults
res_lst[["13G_EXP_ 9.11"]]$"Tumor_only_2_vs_Tumor_only_1_vs_Tumor_only_0" <- res_ts_lst[["13G_EXP_ 9.11"]]
res_lst[["15O_EXP_ 9.11"]]$"Tumor_only_2_vs_Tumor_only_1_vs_Tumor_only_0" <- res_ts_lst[["15O_EXP_ 9.11"]]

#DESeq2 differential expression analysis statistics
res_cline_exper <- lapply(names(dds_cline_exper), function(cline_exper){
  lapply(names(res_lst[[cline_exper]]), function(comparison){
    res_lst[[cline_exper]][[comparison]]%>%
      as.data.frame()%>%
      mutate(comparison = paste(comparison, collapse = "_vs_"),
             data_subset = cline_exper,
             is_sig = ifelse(padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff,TRUE, FALSE),
             is_sig = ifelse(is.na(is_sig), FALSE, is_sig))%>%
      rownames_to_column("barcode_seq")
  }) %>%
    bind_rows() %>%
    arrange( padj, desc(log2FoldChange))
})%>%
  setNames(.,names(dds_cline_exper))

# significant barcodes counts
sig_res_count <- lapply(names(res_cline_exper), function(cline_exper){
  deb <- res_cline_exper[[cline_exper]] %>%
    filter(is_sig) %>% 
    arrange(desc(log2FoldChange), padj) %>%
    pull(barcode_seq) %>%
    unique()
  cts_cline_exper[[cline_exper]] [deb,] 
}) %>%
  setNames(., names(res_cline_exper))

#deseq normalization
deseq_counts_ntd <- lapply(names(dds_cline_exper), function(cline_exper){
  mat <- dds_cline_exper[[cline_exper]]  %>%
    DESeq2::normTransform() %>%
    assay()
  cols <- sample_list[[cline_exper]]
  mat[, cols] 
}) %>%
  setNames(., names(dds_cline_exper))

#significant barcodes normalized counts
sig_res_norm_ntd <- lapply(names(res_cline_exper),function(cline_exper){
  deb <- res_cline_exper[[cline_exper]] %>%
    filter(is_sig) %>% 
    pull(barcode_seq)%>%
    unique()
  deseq_counts_ntd[[cline_exper]][deb, ] 
}) %>%
  setNames(., names(res_cline_exper))

## ---- save_tables
sig_column_header <- glue("Significant (padj < {padj_cutoff} & absolute(log2FoldChange) >= {log2fc_cutoff})")
res_cline_exper %>%
  bind_rows() %>%
  dplyr::select(data_subset, comparison, everything()) %>%
  dplyr::rename(!!sig_column_header := is_sig) %>%
  write_csv(glue("{docs_dir}/deseq2_differential_abundance/deseq2_results_cline_exper_9.11.csv"))

res_cline_exper %>%
  bind_rows() %>%
  dplyr::select(data_subset, comparison, everything()) %>%
  dplyr::count(data_subset,comparison,is_sig) %>%
  dplyr::rename(!!sig_column_header := is_sig) %>%
  write_csv(glue("{docs_dir}/deseq2_differential_abundance/deseq2_results_sig_summary_9.11.csv"))

lapply(deseq_counts_ntd, function(df){
  df %>% as.data.frame %>% rownames_to_column("barcode_seq")
}) %>% 
  openxlsx::write.xlsx(glue("{docs_dir}/counts/cline_exper_blue_aggr_normlog2_9.11.xlsx"))

## Summary of differential expression results
pdf(glue("{plots_dir}/deseq2_results_sig_summary_9.11.pdf"))
res_cline_exper %>%
  bind_rows() %>%
  dplyr::select(data_subset, comparison, everything()) %>%
  dplyr::count(data_subset,comparison,is_sig) %>%
  mutate(t = str_sub(comparison, -1),
         ratio = n/sum(n),
         data_subset = str_replace(data_subset, "_EXP", ""),
         comparison_data = paste0(comparison, "_", data_subset))%>%
  ggplot(aes(is_sig, comparison_data))+
  geom_tile(aes( fill = ratio), show.legend = FALSE)+
  geom_text(aes(label = n))+
  scale_fill_viridis_c()+
  facet_grid(data_subset~., scales = "free_y", space = "free_y")+
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        strip.text= element_text(size = 6),
        panel.spacing = unit(0.1, "cm"))+
  labs(y = "DESeq2 comparisons",
       x = "Is barcode differentially abundant?")

res_cline_exper %>%
  bind_rows() %>%
  dplyr::select(data_subset, comparison, everything()) %>%
  dplyr::count(data_subset,comparison,is_sig) %>%
  mutate(t = paste0("Tumor_only_",str_sub(comparison, -1)),
         data_subset = str_replace(data_subset, "_EXP", ""),
         comparison_data = paste0(comparison, "_", data_subset),
         condition = str_remove(comparison, "_vs.+"))%>%
  group_by(comparison_data) %>% 
  mutate(ratio = n/sum(n))%>%
  ungroup() %>% 
  ggplot(aes(ratio, condition, fill = is_sig, label = n))+
  geom_col()+
  geom_text(position = position_fill(),hjust=1.5)+
  facet_grid(data_subset~t, scales = "free_y", space = "free_y")+
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        strip.text= element_text(size = 6),
        panel.spacing = unit(0.1, "cm"))+
  labs(y = "Condition",
       x = "Group ratio")

dev.off()


##########
# Explore results
# Explore differentially abundant barcodes
##Ordered heatmaps for barcodes from each comparison showing all the samples
pdf(glue("{plots_dir}/dab_normalized_compraisons_split_padj{padj_cutoff}_log2fc{log2fc_cutoff}_all_samples_9.11.pdf"))
lapply(names(res_cline_exper), function(res_name){
  #select differentially abundant barcodes
  comparisons_lst <- res_cline_exper[[res_name]] %>%
    filter(is_sig) %>%
    dplyr::select(barcode_seq, comparison) %$%
    split(barcode_seq, comparison)
  cols <- colnames(sig_res_norm_ntd[[res_name]])
  hm_mat <- t(apply(sig_res_norm_ntd[[res_name]],1,scale))
  barcodes_n <- nrow(cts_cline_exper[[res_name]])
  colnames(hm_mat) <- cols
  #plot all samples
  lapply(names(comparisons_lst), function(comparison_x){
    comparison_dab <- comparisons_lst[[comparison_x]]
    ComplexHeatmap::Heatmap(hm_mat[comparison_dab, ], cluster_columns=FALSE,row_title_gp = gpar( fontsize = 3),
                            column_names_gp = gpar(fontsize = 8),
                            show_row_names=FALSE,
                            show_column_names=TRUE,
                            column_title = glue("{res_name} - {comparison_x}\npadj < {padj_cutoff} & absolute(log2FoldChange) >= {log2fc_cutoff} ({length(comparison_dab)}/{barcodes_n})"))
  })
})
dev.off()

##Clustered heatmaps for barcodes from each comparison showing all the samples
pdf(glue("{plots_dir}/dab_normalized_compraisons_split_padj{padj_cutoff}_log2fc{log2fc_cutoff}_all_samples_9.11_clustered.pdf"))
lapply(names(res_cline_exper), function(res_name){
  #select differentially abundant barcodes
  comparisons_lst <- res_cline_exper[[res_name]] %>%
    filter(is_sig) %>%
    dplyr::select(barcode_seq, comparison) %$%
    split(barcode_seq, comparison)
  cols <- colnames(sig_res_norm_ntd[[res_name]])
  hm_mat <- t(apply(sig_res_norm_ntd[[res_name]],1,scale))
  barcodes_n <- nrow(cts_cline_exper[[res_name]])
  colnames(hm_mat) <- cols
  #plot all samples
  lapply(names(comparisons_lst), function(comparison_x){
    comparison_dab <- comparisons_lst[[comparison_x]]
    message(length(comparison_dab))
    cols <- colnames(hm_mat)
    column_ha <- HeatmapAnnotation(foo1 = str_extract(cols, "Tumor_only|Tumor\\+NKcells|Tumor\\+PerforinKONKcells|Tumor\\+IFNgKONKcells"),
                                   bar1 =  str_extract(cols, "timepoint."),
                                   col = list(foo1 = c("Tumor_only" = "#999999", "Tumor+IFNgKONKcells" = "#984ea3", "Tumor+NKcells" = "#e41a1c", "Tumor+PerforinKONKcells" = "#ff7f00"),
                                              bar1 = c("timepoint0" = "#e5f5f9", "timepoint1" = "#41ae76", "timepoint2" = "#00441b")))
    ComplexHeatmap::Heatmap(hm_mat[comparison_dab, ],
                            cluster_columns=TRUE,row_title_gp = gpar( fontsize = 3),
                            column_names_gp = gpar(fontsize = 8),
                            show_row_names=FALSE,
                            show_column_names=TRUE,
                            top_annotation = column_ha, 
                            column_title = glue("{res_name} - {comparison_x}\npadj < {padj_cutoff} & absolute(log2FoldChange) >= {log2fc_cutoff} ({length(comparison_dab)}/{barcodes_n})"))
  })
})
dev.off()

# Split output into lists
##split differential abundance results by comparisons and cell lines
comp_lst_all <- res_cline_exper %>%
  bind_rows(.id = "cline") %>%
  split(.,.$comparison) %>%
  lapply(function(x) split(x,x$cline))

##list of samples at reference timepoint and timepoint 2
t_cols_lst <- lapply(names(comp_lst_all), function(comp){
  t_comps <- str_split(comp, "_vs_")[[1]]
  t_ref <- t_comps[length(t_comps)]
  t_2 <- t_comps[1]
  cline_exper_lst <- comp_lst_all[[comp]] 
  lapply(names(cline_exper_lst), function(cline){
    #reference samples
    tref_cols <- coldata %>%
      mutate(cline_exper = paste0(cell_line,"_", experiment)) %>%
      filter(cline_exper== cline & sample_bio_t == t_ref) %>%
      rownames()
    #t2 samples
    t2_cols <- coldata %>%
      mutate(cline_exper = paste0(cell_line,"_", experiment)) %>%
      filter(cline_exper== cline & sample_bio_t == t_2) %>%
      rownames()
    l <- list(tref_cols=tref_cols, t2_cols=t2_cols)
    l
  }) %>%
    setNames(., names(cline_exper_lst))
})%>%
  setNames(., names(comp_lst_all))

#barcode variability
barcode_var <- function(x){
  x <- x+1
  log2(max(x)/min(x))
}

#Define barcode groups
cline_exper_barcode_groups <- lapply(names(t_cols_lst), function(comp){
  comp_lst_x <- t_cols_lst[[comp]]
  
  lapply(names(comp_lst_x), function(res_name){
    #normalized matrix
    hm_norm <- deseq_counts_ntd[[res_name]]
    #count matrix
    hm_count <- cts_cline_exper[[res_name]]
    tref_cols <- comp_lst_x[[res_name]][["tref_cols"]]
    t2_cols <- comp_lst_x[[res_name]][["t2_cols"]]
    cols <- c(tref_cols, t2_cols)
    #scale the normalized sub-matrix of samples of interest
    hm_norm_scaled <- t(apply(hm_norm[,cols],1,scale))
    colnames(hm_norm_scaled) <- cols
    hm_norm_scaled <- hm_norm_scaled[,cols]
    #Deseq statistics of the comparison of interest 
    comparison_dab <- comp_lst_all[[comp]][[res_name]] %>%
      arrange(padj) %>%
      mutate(dir = ifelse(log2FoldChange>=1,"increasing",ifelse(log2FoldChange<=-1,"decreasing",NA))) #define increasing and decreasing barcode groups
    #variability measures at tref and t2
    tref_rangelog2fc <- apply(hm_norm[comparison_dab$barcode_seq, tref_cols], 1, barcode_var)
    t2_rangelog2fc <- apply(hm_norm[comparison_dab$barcode_seq, t2_cols], 1, barcode_var)
    #add variability measures to the camparison table for easy ordering
    comparison_dab <- cbind(comparison_dab,
                            tref_rangelog2fc,
                            t2_rangelog2fc)
    #define timepoints cutoffs
    tref_cutoff <- 0.5
    t2_cutoff <- 1
    
    comparison_dab %>%
      as.data.frame %>%
      mutate(barcode_group = case_when(
        is_sig & dir == "increasing" ~ "primary_resistant",
        is_sig & dir == "decreasing" ~ "eliminated",
        tref_rangelog2fc <= tref_cutoff & t2_rangelog2fc >= t2_cutoff ~ "immunoedited",
        tref_rangelog2fc <= tref_cutoff & !is.na(dir) ~ "others",
        tref_rangelog2fc <= tref_cutoff & t2_rangelog2fc <= t2_cutoff ~ "static",
        TRUE ~ "others")) 
    
  }) %>%
    setNames(., names(comp_lst_x))
  
}) %>%
  setNames(., names(t_cols_lst))

#
define_immunoedit_groups <- function(t2_df_sorted, t2_n_sorted){
  ##cluster barcodes after sorting
  sort_cor_dist <- 1-cor(t2_df_sorted %>% t() %>% scale())
  hc_2 <- hclust(as.dist(sort_cor_dist), method = "ward.D2")
  groups_list_sorted <- cutree(hc_2,t2_n_sorted) %>% split(names(.),.)
  
  ##assign group identity to the clusters based on the scaled abundance of middle well
  t2_group_names <- c("secondary_resistant", "negatively_edited")
  t2_group_order <- sapply(groups_list_sorted,function(grp){
    df <- t2_df_sorted %>% t() %>% scale() %>% t()
    mean(df[grp,2])
  }) %>% order
  names(groups_list_sorted) <- t2_group_names[t2_group_order]
  group_df <- groups_list_sorted %>%
    stack() %>%
    as.data.frame() %>%
    dplyr::rename("barcode_seq" = "values",
                  "t2_group_sorted" = "ind")
  return(group_df)
}

#Create R object with all barcodes information
barcode_info <- lapply(names(t_cols_lst), function(comparison_x){
  comp_lst_x <- t_cols_lst[[comparison_x]]
  lapply(names(comp_lst_x), function(res_name){
    #normalized matrix
    hm_norm <- deseq_counts_ntd[[res_name]]
    #count matrix
    hm_count <- cts_cline_exper[[res_name]]
    #column names 
    tref_cols <- comp_lst_x[[res_name]][["tref_cols"]]
    t2_cols <- comp_lst_x[[res_name]][["t2_cols"]]
    
    #DESeq statistics of the comparison of interest 
    group_3_immun <- cline_exper_barcode_groups[[comparison_x]][[res_name]]
    
    #Immunoedited barcode clustering to define positively and negative edited barcodes
    ##Sort barcode abundances at t2 samples so that wells with highest abundances end-up in the same column
    immunoedit_barcodes <- filter(group_3_immun, barcode_group == "immunoedited")$barcode_seq
    t2_df <- hm_norm[immunoedit_barcodes,t2_cols]
    t2_n_sorted <- 2 #the "secondary_resistant" and "negatively_edited" subgroups
    t2_df_sorted <- apply(t2_df,1,sort) %>% t()
    colnames(t2_df_sorted) <- colnames(t2_df)
    group_df <- define_immunoedit_groups(t2_df_sorted, t2_n_sorted)
    
    #normalize values by the average t0 abundance
    t2_group_norm <- hm_norm[group_df$barcode_seq, c(tref_cols,t2_cols)]
    t2_group_norm_t0_scaled <- t2_group_norm %>%
      apply(1, function(x){
        x/mean(x[1:3])
      }) %>%
      t()
    #variability at t2 is increasing or not
    group_df$t2_max_well_fc <- apply(t2_group_norm_t0_scaled[,t2_cols], 1, max) #maximum value at t2 after dividing by t0 average
    group_df$t2_rangelog2fc <- group_3_immun[group_df$barcode_seq,"t2_rangelog2fc"]
    
    #add defined barcode clusters assignment
    group_3_immun[group_df$barcode_seq,"barcode_group"] <- as.character(group_df$t2_group_sorted)
    group_3_immun[group_df$barcode_seq,"t2_max_well_fc"] <- group_df$t2_max_well_fc
    group_3_immun <- group_3_immun %>%
      arrange(barcode_group) %>%
      mutate(barcode_group = as.factor(barcode_group),
             barcode_group = forcats::fct_relevel(barcode_group, "primary_resistant", "eliminated"),
             barcode_group = forcats::fct_relevel(barcode_group, "static", "others", after = Inf),
             barcode_group_2 = ifelse(str_detect(barcode_group, "edited"), "immunoedited", as.character(barcode_group)),
             barcode_group_2 = forcats::fct_relevel(barcode_group_2, "primary_resistant", "eliminated", "static", "immunoedited", "others"))%>%
      arrange(barcode_group)
    
    #save all the objects
    list(cols = colnames(hm_norm),
         hm_count = hm_count, 
         hm_norm = hm_norm,
         hc_2 = hc_2,
         t2_names = t2_group_names[t2_group_order],
         group_df = group_df,
         t2_df = t2_df[group_df$barcode_seq,],
         t2_df_sorted = t2_df_sorted[group_df$barcode_seq,],
         t2_group_norm = t2_group_norm[group_df$barcode_seq,],
         t2_group_norm_t0_scaled = t2_group_norm_t0_scaled[group_df$barcode_seq,],
         group_3_immun = group_3_immun
    )
  }) %>% setNames(names(comp_lst_x)) 
}) %>% setNames(names(t_cols_lst)) 

#plot immunoedited barcodes variability
scatter_immunoedited_var <- function(barcode_groups_obj){
  comp_names <- names(barcode_groups_obj)
  lapply(comp_names, function(comp){
    cline_names <- barcode_groups_obj[[comp]] %>% names
    lapply(cline_names, function(cline){
      comp_cline_df <- barcode_groups_obj[[comp]][[cline]]$group_3_immun
      comp_cline_df %>% 
        filter(barcode_group %in% c("others", "secondary_resistant", "static")) %>% 
        ggplot(aes(tref_rangelog2fc,t2_rangelog2fc, color = barcode_group))+
        geom_point()+
        labs(title = glue::glue("{comp} {cline}"))
    })
  })
}

pdf(glue("{plots_dir}/scatter_immunoedited_var_compraisons_cline.pdf"))
scatter_immunoedited_var(barcode_info)
dev.off()

#Combine output statistics and defined barcode groups across comparisons of each cell line
barcode_info_flat <- barcode_info %>% unlist(recursive=FALSE) 
## shorten comparisons' names
short_names <- barcode_info_flat %>% names() %>% str_replace_all("umor|erforin|cells_", "") %>% str_replace_all("_only_","only")
## save barcode raw counts,comparisons statistics, and assigned group
lapply(barcode_info_flat, function(x){
  x[["hm_count"]] %>% 
    rownames_to_column("barcode_seq") %>% 
    left_join(x[["group_3_immun"]])
}) %>%
  setNames(., short_names) %>% 
  openxlsx::write.xlsx(glue("{docs_dir}/counts/barcodes_counts_comparisonstats_groups.xlsx"))

## Per comparison, save barcode normalized counts,comparisons statistics, and assigned group
lapply(barcode_info_flat, function(x){
  x[["hm_norm"]] %>% 
    as.data.frame() %>% 
    rownames_to_column("barcode_seq") %>% 
    left_join(x[["group_3_immun"]])
}) %>%
  setNames(., short_names) %>% 
  openxlsx::write.xlsx(glue("{docs_dir}/counts/barcodes_norm_comparisonstats_groups.xlsx"))

## Per comparison comparisons statistics, and assigned group
lapply(barcode_info_flat, function(x){
  x[["group_3_immun"]]
}) %>%
  setNames(., short_names) %>% 
  openxlsx::write.xlsx(glue("{docs_dir}/counts/barcodes_comparisons_stats_groups.xlsx"))

## Per cell line group, save barcode normalized counts and assigned group in each comparison
cline_comps_lst <- barcode_info_flat %>% names %>% split(.,str_extract(., "(13G|15O|15M|13H).+"))

groups_lst <- lapply(cline_comps_lst,function(cline_comp){
  lapply(cline_comp, function(comp){
    col_name <- glue("{comp}")
    comp_lst <- barcode_info_flat[[comp]]
    comp_lst[["group_3_immun"]] %>% 
      dplyr::select(c("barcode_seq","barcode_group_2")) %>% 
      dplyr::rename(!!col_name := barcode_group_2)
  }) %>% 
    purrr::reduce(full_join)
}) 

lapply(names(groups_lst), function(cline){
  deseq_counts_ntd[[cline]]%>% 
    as.data.frame() %>% 
    rownames_to_column("barcode_seq") %>% 
    full_join(groups_lst[[cline]]) 
}) %>% 
  setNames(., names(groups_lst)) %>% 
  openxlsx::write.xlsx(glue("{docs_dir}/counts/cline_norm_groups.xlsx"))

lapply(names(groups_lst), function(cline){
  cts_cline_exper[[cline]] %>% 
    rownames_to_column("barcode_seq") %>% 
    full_join(groups_lst[[cline]])%>% 
    as.data.frame()
}) %>% 
  setNames(., names(groups_lst)) %>% 
  openxlsx::write.xlsx(glue("{docs_dir}/counts/cline_counts_groups.xlsx"))

##comparisons statistics and barocode groups
res_name_samples <- barcode_info_flat %>% names %>% str_replace("^.+only_[0-9]\\.","")
res_name_info_lst <- split(barcode_info_flat, res_name_samples)
stats_lst <- lapply(names(res_name_info_lst), function(res_name){
  lst <- res_name_info_lst[[res_name]]
  names(lst) <- names(lst) %>% str_replace(paste0(".",res_name), "")  
  barcode_groups_stats <- lapply(names(lst), function(comparison_x){
    #barcodes information
    comparison_info_x <- lst[[comparison_x]]
    #extract experiment name
    exper <- str_split(res_name,"_",simplify=TRUE)[,3] %>% unique %>% paste(collapse = "")
    exper <- ifelse(exper == "9.8", "9.8", "9.67")
    
    #select barcode group info of interest
    barcode_groups <- comparison_info_x$group_3_immun %>%
      dplyr::select(barcode_seq, tref_rangelog2fc, t2_rangelog2fc, barcode_group)
    barcode_groups %>%  rename_at(vars(-barcode_seq), ~paste0(.x,"_", comparison_x)) 
  }) %>%
    purrr::reduce(by = "barcode_seq", full_join)
  
  #select deseq statistics  
  barcode_deseq_stats <-res_cline_exper[[res_name]] %>% 
    dplyr::select(barcode_seq,log2FoldChange, padj, comparison) %>%
    pivot_wider(names_from = comparison,
                values_from = c(log2FoldChange, padj))
  
  #merge groups and deseq statistics 
  left_join(barcode_groups_stats,barcode_deseq_stats) %>%
    dplyr::select(barcode_seq, contains("barcode_group"), everything())
}) %>%
  setNames(names(res_name_info_lst))
#Save stats and barcodes groups of each cell lines 
openxlsx::write.xlsx(stats_lst,
                     glue("{docs_dir}/counts/barcode_groups_stats.xlsx"))

#heatmap of barcodes numbers
pdf(glue("{plots_dir}/number_and_percentage_of_barcodes_per_group_9.11_ref.pdf"), width = 10)
map_depth(barcode_info,2,`[[`,"group_3_immun") %>%
  lapply(bind_rows, .id = "cline") %>%
  bind_rows(.id = "comp") %>%
  dplyr::count(comp, cline, barcode_group_2) %>%
  group_by(cline, comp) %>%
  mutate(percent = n/sum(n),
         n2 = ifelse(barcode_group_2 != "secondary_resistant", 0, n)) %>%
  arrange(cline,desc(n2)) %>% 
  ungroup() %>% 
  mutate(comp = factor(comp, levels = rev(unique(comp)))) %>% 
  ggplot( aes(barcode_group_2, comp ,fill = percent, label = n))+
  geom_tile()+
  geom_text(color = "white")+
  scale_fill_viridis_c()+
  facet_grid(cline~., scales = "free")+
  labs(x = "group", y = "", fill = "barcode percentage\nper group", title  = glue::glue("Number of barcodes: {comparison_x}"))+
  theme(title = element_text(size = 8),
        axis.text.x = element_text(size = 10, angle = 90)) 
dev.off()

#reads per group
groups_sample_reads <- lapply(names(barcode_info), function(comparison_x){
  comparison_info_x <- barcode_info[[comparison_x]]
  lapply(names(comparison_info_x), function(res_name){
    #barcode info
    cline_info <- comparison_info_x[[res_name]]
    #split barcodes by group
    groups_lst <- cline_info[["group_3_immun"]] %$%
      split(barcode_seq, barcode_group_2)
    #sum number of barcode per group
    groups_sample_reads <- lapply(groups_lst, function(group_x){
      colSums(cline_info$hm_count[group_x,]) %>%
        enframe() 
    }) %>%
      bind_rows(.id = "barcode_group") %>%
      mutate(barcode_group = factor(barcode_group,
                                    levels = c("primary_resistant", "eliminated", "static","immunoedited","others", "secondary_resistant")))
  })%>%
    setNames(names(comparison_info_x))
})%>%
  setNames(names(barcode_info))

groups_sample_reads %>%
  lapply(bind_rows, .id = "cline") %>%
  bind_rows( .id = "comparison") %>%
  dplyr::rename(n_reads = "value") %>%
  write_csv(glue("{docs_dir}/counts/cline_number_of_reads_per_group_per_comparison_9.11.csv"))

#heatmap of reads numbers
pdf(glue("{plots_dir}/heatmaps_number_of_reads_per_sample_9.11_ref.pdf"))
##separate plots for each cell-line experiment 
lapply(names(groups_sample_reads), function(comparison_x){
  lapply(names(barcode_info[[comparison_x]]), function(res_name){
    groups_sample_reads[[comparison_x]][[res_name]] %>%
      group_by(name) %>%
      mutate(percent = value/sum(value),
             cell_line = res_name) %>%
      ungroup()
  }) %>%
    lapply(function(df){
      df %>%
        ggplot(aes(barcode_group, name ,fill = percent, label = value))+
        geom_tile()+
        geom_text(color = "white", size = 3)+
        scale_fill_viridis_c()+
        labs(x = "group", y = "", fill = "group\npercentage", title = glue::glue("{unique(df$cell_line)} {comparison_x}: number of reads"))+
        theme(axis.text.x = element_text(size = 10, angle = 90),
              axis.text.y = element_text(size = 6),
              title = element_text(size = 8))
    })
})
dev.off()

#percentage of reads per group per sample
pdf(glue("{plots_dir}/barplots_number_of_reads_per_sample_9.11_ref.pdf"))
lapply(names(groups_sample_reads), function(comparison_x){
  lapply(names(barcode_info[[comparison_x]]), function(res_name){
    groups_sample_reads[[comparison_x]][[res_name]] %>%
      group_by(name) %>%
      mutate(percent = value/sum(value)) %>%
      ungroup()%>%
      ggplot(aes(percent, name, fill = barcode_group))+
      geom_col(position = position_fill())+
      scale_fill_brewer(palette ="Dark2", direction = 1)+
      theme(axis.text.y = element_text(size = 4))+
      labs(x = "reads percentage", title = glue::glue("{res_name} {comparison_x}: number of reads"))+
      theme(title = element_text(size = 8))
  })
}) 
dev.off()

#histogram of distribution of reads per well per group
pdf(glue("{plots_dir}/barcodes_count_distribution_group_wells_9.11_ref.pdf"),width=12,height = 5)
barcodes_count_group_wells <- lapply(names(groups_sample_reads), function(comparison_x){
  lapply(names(barcode_info[[comparison_x]]), function(res_name){
    barcode_info[[comparison_x]][[res_name]]$hm_count %>%
      rownames_to_column("barcode_seq") %>%
      left_join(barcode_info[[comparison_x]][[res_name]]$group_3_immun %>% dplyr::select(barcode_seq, barcode_group, barcode_group_2)) %>%
      dplyr::select(barcode_seq,barcode_group_2,everything(), -barcode_group) %>%
      pivot_longer(-c(barcode_seq,barcode_group_2)) 
  }) %>%
    setNames(.,names(barcode_info[[comparison_x]]))
})%>%
  setNames(.,names(groups_sample_reads))

lapply(names(barcodes_count_group_wells), function(comparison_x){
  lapply(names(barcode_info[[comparison_x]]), function(res_name){
    barcodes_count_group_wells[[comparison_x]][[res_name]] %>%
      filter(value != 0)%>%
      mutate(name = str_extract(name, "T.+")) %>% 
      ggplot(aes(value))+
      geom_histogram()+
      scale_x_log10()+
      facet_grid(barcode_group_2~name, scales = "free_y")+
      labs(title = res_name, x = "Number of reads per barcode", y = "Frequency")+
      theme(axis.text.x = element_text(size = 7),
            axis.text.y = element_text(size = 7),
            strip.text.x = element_text(size = 5),
            strip.text.y = element_text(size = 6))
  })
})
dev.off()

#heatmaps
pdf(glue("{plots_dir}/barcode_info_heatmaps_9.11_ref.pdf"),width = 12)
walk(names(barcode_info), function(comparison_x){
  comparison_info_x <- barcode_info[[comparison_x]]
  lapply(names(comparison_info_x), function(cline){
    #get all the information about the cell line and comparison of interest
    info <- comparison_info_x[[cline]]
    #get ordered column names
    cols_lst <- t_cols_lst[[comparison_x]][[cline]] #%>% unlist()
    cond_comp <- cols_lst$t2_cols %>% str_remove("_timepoint.+") %>% unique
    cols_all <- info$cols
    cols_comp <- c(cols_lst$tref_cols, str_subset(cols_all, stringr::fixed(cond_comp))) %>% unique()
    #heatmap of the sorted t2 wells
    hm_cor_df <- info$t2_df_sorted %>% t() %>% scale() %>% t()
    hm_cor <- ComplexHeatmap::Heatmap(hm_cor_df,
                                      #cluster_rows = info$hc_2,
                                      row_split = info$group_df$t2_group_sorted %>% as.character() %>% str_replace("_[a-z-]+$",""),
                                      row_title_rot = 0,
                                      cluster_columns=FALSE,
                                      show_row_names=FALSE,
                                      name = "scaled log2(counts+1)",
                                      row_dend_reorder = TRUE,
                                      cluster_row_slices = FALSE,
                                      column_title = glue("t2_non-constant barcodes at t2\n {comparison_x} {cline}"),
                                      width = ncol(info$t2_df_sorted)*unit(6, "mm")) 
    hm_cor <- ComplexHeatmap::draw(hm_cor)
    #label all the barcodes into the different defined groups
    hm_row_split_labels <- info$group_3_immun %>%
      mutate(barcode_group_2_a = factor(barcode_group_2,
                                        levels = c("primary_resistant", "eliminated", "static", "immunoedited", "others", "secondary_resistant")))%>%
      group_by(barcode_group_2) %>%
      mutate(barcode_group_2 = glue("{barcode_group_2} (n={n()})")) %>%
      ungroup()%>%
      arrange(barcode_group_2_a)%>%
      dplyr::select(barcode_seq,barcode_group_2) %>%
      mutate(barcode_group_2 = factor(barcode_group_2,
                                      levels = unique(barcode_group_2)))
    #row annotations
    t1_counts <- apply(info$hm_count[hm_row_split_labels$barcode_seq,cols_all][,1:3], 1, mean)
    t1_counts_3 <- apply(info$hm_count[filter(hm_row_split_labels,!str_detect(barcode_group_2,"others|immunoedited"))$barcode_seq,cols_all][,1:3], 1, mean)
    
    row_annot <- ComplexHeatmap::HeatmapAnnotation(which = "row",
                                                   `Avg count t0\n(log10+1)` =  ComplexHeatmap::anno_barplot(log10(t1_counts+1),axis_param = list(side = "top")),
                                                   show_annotation_name = TRUE,
                                                   annotation_name_gp = gpar(fontsize = 6 , fontface = "bold"))
    row_annot_3 <- ComplexHeatmap::HeatmapAnnotation(which = "row",
                                                     `Avg count t0\n(log10+1)` =  ComplexHeatmap::anno_barplot(log10(t1_counts_3+1),axis_param = list(side = "top")),
                                                     show_annotation_name = TRUE,
                                                     annotation_name_gp = gpar(fontsize = 6 , fontface = "bold"))
    #column annotation
    ##all samples annotation
    hm_sh_div <- shannon_diversity_lst[[cline]][cols_all]
    column_annot <- ComplexHeatmap::HeatmapAnnotation(`Shannon diversity` =  ComplexHeatmap::anno_lines(hm_sh_div, add_points = TRUE, ylim = c(0,max(hm_sh_div)+(diff(range(hm_sh_div))*0.2))),
                                                      condition = setNames(coldata$sample_bio,rownames(coldata))[cols_all],
                                                      time_point = setNames(coldata$timepoint,rownames(coldata))[cols_all],
                                                      show_legend = TRUE,
                                                      simple_anno_size = unit(0.35, "cm"),
                                                      show_annotation_name = TRUE,
                                                      col = list(condition = c("Tumor_only" = "#35978f","Tumor+Perforin KO NK cells" = "#dfc27d", "Tumor+IFNg KO NK cells" = "#bf812d", "Tumor+NK cells" = "#8c510a"),
                                                                 time_point = c("timepoint_0" = "#d9d9d9", "timepoint_1" = "#969696", "timepoint_2" = "#525252")),
                                                      annotation_name_gp = gpar(fontsize = 8, fontface = "bold"))
    #all groups
    hm1_mat <- info$hm_norm[hm_row_split_labels$barcode_seq,cols_all]%>%t%>%scale()%>%t()
    hm1 <- ComplexHeatmap::Heatmap(hm1_mat,
                                   name = "scaled log2(counts+1)",
                                   cluster_columns=FALSE,
                                   cluster_rows=FALSE,
                                   row_split = hm_row_split_labels$barcode_group_2,
                                   row_title_rot = 0,
                                   row_title_gp = gpar( fontsize = 8),
                                   show_row_names=FALSE,
                                   show_column_names=TRUE,
                                   top_annotation = column_annot,
                                   left_annotation = row_annot,
                                   column_names_gp = gpar(fontsize = 8),
                                   column_title = glue("{comparison_x} {cline}"),
                                   width = length(cols_all)*unit(6, "mm")) 
    ComplexHeatmap::draw(hm1,
                         merge_legend = TRUE)
    #exclude others and immunoedited
    hm3_mat <- info$hm_norm[filter(hm_row_split_labels,!str_detect(barcode_group_2,"others|immunoedited"))$barcode_seq,cols_all]%>%t%>%scale()%>%t()
    hm3 <- ComplexHeatmap::Heatmap(hm3_mat,
                                   name = "scaled log2(counts+1)",
                                   cluster_columns=FALSE,
                                   cluster_rows=FALSE,
                                   row_split = filter(hm_row_split_labels,!str_detect(barcode_group_2,"others|immunoedited"))$barcode_group_2,
                                   row_title_rot = 0,
                                   row_title_gp = gpar( fontsize = 8),
                                   show_row_names=FALSE,
                                   show_column_names=TRUE,
                                   top_annotation = column_annot,
                                   left_annotation = row_annot_3,
                                   column_names_gp = gpar(fontsize = 8),
                                   column_title = glue("{comparison_x} {cline}"),
                                   width = length(cols_all)*unit(6, "mm")) 
    ComplexHeatmap::draw(hm3,
                         merge_legend = TRUE)
    
    
    #comparison samples
    hm_sh_div_comp <- shannon_diversity_lst[[cline]][cols_comp]
    column_annot_comp <- ComplexHeatmap::HeatmapAnnotation(`Shannon diversity` =  ComplexHeatmap::anno_lines(hm_sh_div_comp, add_points = TRUE, ylim = c(0,max(hm_sh_div_comp)+(diff(range(hm_sh_div_comp))*0.2))),
                                                           condition = setNames(coldata$sample_bio,rownames(coldata))[cols_comp],
                                                           time_point = setNames(coldata$timepoint,rownames(coldata))[cols_comp],
                                                           show_legend = TRUE,
                                                           simple_anno_size = unit(0.35, "cm"),
                                                           show_annotation_name = TRUE,
                                                           col = list(condition = c("Tumor_only" = "#35978f","Tumor+Perforin KO NK cells" = "#dfc27d", "Tumor+IFNg KO NK cells" = "#bf812d", "Tumor+NK cells" = "#8c510a"),
                                                                      time_point = c("timepoint_0" = "#d9d9d9", "timepoint_1" = "#969696", "timepoint_2" = "#525252")),
                                                           annotation_name_gp = gpar(fontsize = 8, fontface = "bold"))
    
    ##all groups
    hm1_mat_comp <- info$hm_norm[hm_row_split_labels$barcode_seq,cols_comp]%>%t%>%scale()%>%t()
    hm1_comp <- ComplexHeatmap::Heatmap(hm1_mat_comp,
                                        name = "scaled log2(counts+1)",
                                        cluster_columns=FALSE,
                                        cluster_rows=FALSE,
                                        row_split = hm_row_split_labels$barcode_group_2,
                                        row_title_rot = 0,
                                        row_title_gp = gpar( fontsize = 8),
                                        show_row_names=FALSE,
                                        show_column_names=TRUE,
                                        top_annotation = column_annot_comp,
                                        left_annotation = row_annot,
                                        column_names_gp = gpar(fontsize = 8),
                                        column_title = glue("{comparison_x} {cline}"),
                                        width = length(cols_comp)*unit(6, "mm")) 
    ComplexHeatmap::draw(hm1_comp,
                         merge_legend = TRUE)
    ##exclude others and immunoedited
    hm3_mat_comp <- info$hm_norm[filter(hm_row_split_labels,!str_detect(barcode_group_2,"others|immunoedited"))$barcode_seq,cols_comp]%>%t%>%scale()%>%t()
    hm3_comp <- ComplexHeatmap::Heatmap(hm3_mat_comp,
                                        name = "scaled log2(counts+1)",
                                        cluster_columns=FALSE,
                                        cluster_rows=FALSE,
                                        row_split = filter(hm_row_split_labels,!str_detect(barcode_group_2,"others|immunoedited"))$barcode_group_2,
                                        row_title_rot = 0,
                                        row_title_gp = gpar( fontsize = 8),
                                        show_row_names=FALSE,
                                        show_column_names=TRUE,
                                        top_annotation = column_annot_comp,
                                        left_annotation = row_annot_3,
                                        column_names_gp = gpar(fontsize = 8),
                                        column_title = glue("{comparison_x} {cline}"),
                                        width = length(cols_comp)*unit(6, "mm")) 
    ComplexHeatmap::draw(hm3_comp,
                         merge_legend = TRUE)
  })
  
})
dev.off()


pdf(glue("{plots_dir}/barcode_info_bubble_9.11_ref_unsorted.pdf"), height = 8, width = 8)
walk(names(barcode_info), function(comparison_x){
  comparison_info_x <- barcode_info[[comparison_x]]
  lapply(names(comparison_info_x), function(res_name){
    info <- comparison_info_x[[res_name]]
    t_ref <- str_sub(comparison_x, nchar(comparison_x), nchar(comparison_x))
    # Reorder replicates at t2 based on the normalized value 
    timepoint_ref <- glue("Tumor_only_timepoint{t_ref}")
    t2_var_line_df <- info$t2_group_norm %>%
      as.data.frame() %>%
      rownames_to_column("barcode_seq")%>%
      left_join(info$group_3_immun) %>%
      pivot_longer(-c(barcode_seq,t2_max_well_fc,cline:barcode_group_2))
    # Add a column with the size of immunoedited subgroups (positive, negative) to the table
    ggdf <- t2_var_line_df %>%
      group_by(barcode_group) %>%
      mutate(barcode_group = glue("{barcode_group} (n = {length(unique(barcode_seq))})")) %>%
      ungroup
    # Total number of immunoedited barcodes 
    barcodes <- length(unique(ggdf$barcode_seq))
    
    #
    ggplt <- ggdf %>%
      mutate(name = str_replace_all(name, "3.15_EXP_9.678_(A|D)_|timepoint",""),
             name = fct_inorder(name),
             barcode_group = as.character(barcode_group)) %>%
      split(., .$barcode_group) %>%
      lapply(function(x){
        ggplt <- x %>% ggplot(aes(name, barcode_seq,
                                  size = value,
                                  color = value))+
          geom_point()+
          scale_radius(limits = c(0, NA), range = c(0, 5))+
          scale_color_viridis_c()+
          labs(y = "Barcode",
               title = glue::glue("{res_name} {comparison_x}"),
               subtitle = glue::glue("{unique(x$barcode_group)}"),
               size = "log2(counts+1)",
               color = "log2(counts+1)")+
          theme(axis.text.x = element_text(angle=90, hjust = 1),
                axis.text.y = element_text(size = 4),
                legend.position = "right",
                panel.grid.major.y = element_line(colour = "grey95"))
        ggplt%>% print()
      })
    #ggplt %>% print()
  })
})
dev.off()


pdf(glue("{plots_dir}/barcode_info_bubble_9.11_ref_unsorted_notext.pdf"), height = 3.5, width = 2.5)
walk(names(barcode_info), function(comparison_x){
  comparison_info_x <- barcode_info[[comparison_x]]
  lapply(names(comparison_info_x), function(res_name){
    info <- comparison_info_x[[res_name]]
    t_ref <- str_sub(comparison_x, nchar(comparison_x), nchar(comparison_x))
    # Reorder replicates at t2 based on the normalized value 
    timepoint_ref <- glue("Tumor_only_timepoint{t_ref}")
    t2_var_line_df <- info$t2_group_norm %>%
      as.data.frame() %>%
      rownames_to_column("barcode_seq")%>%
      left_join(info$group_3_immun) %>%
      pivot_longer(-c(barcode_seq,t2_max_well_fc,cline:barcode_group_2))
    # Add a column with the size of immunoedited subgroups (positive, negative) to the table
    ggdf <- t2_var_line_df %>%
      group_by(barcode_group) %>%
      mutate(barcode_group = glue("{barcode_group} (n = {length(unique(barcode_seq))})")) %>%
      ungroup
    # Total number of immunoedited barcodes 
    barcodes <- length(unique(ggdf$barcode_seq))
    
    #
    ggplt <- ggdf %>%
      mutate(name = str_replace_all(name, "3.15_EXP_9.678_(A|D)_|timepoint",""),
             name = fct_inorder(name),
             barcode_group = as.character(barcode_group)) %>%
      split(., .$barcode_group) %>%
      lapply(function(x){
        ggplt <- x %>% ggplot(aes(name, barcode_seq,
                                  size = value,
                                  color = value))+
          geom_point()+
          scale_radius(limits = c(0, NA), range = c(0, 5))+
          scale_color_viridis_c()+
          labs(y = "Barcode",
               title = glue::glue("{res_name} {comparison_x}"),
               subtitle = glue::glue("{unique(x$barcode_group)}"),
               size = "log2(counts+1)",
               color = "log2(counts+1)")+
          theme(axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                legend.position = "none",
                axis.title = element_blank(),
                panel.grid.major.y = element_line(colour = "grey95"),
                plot.title = element_text(size = 4),
                plot.subtitle = element_text(size = 4))
        
        ggplt%>% print()
      })
    #ggplt %>% print()
  })
})
dev.off()

# common barcodes
## list of barcodes groups per comparison-cell line
barcodeseq_group_lst <- lapply(barcode_info, function(comp_lst){
  group_lst <- lapply(comp_lst,`[[`, "group_3_immun")
  lapply(group_lst, function(cline){
    cline %$% 
      split(barcode_seq, barcode_group)
  }) %>% unlist(recursive=FALSE)
}) %>% unlist(recursive=FALSE)
## split barcodes by cell line
subgroup_names <- names(barcodeseq_group_lst) %>%
  str_split("\\.", simplify= TRUE) 
celline_group_lst <- split(barcodeseq_group_lst, apply(subgroup_names[,2:4], 1, paste, collapse = "."))
celline_group_lst <- celline_group_lst[lengths(celline_group_lst)>1]

## list  of matrices of the number of common barcodes
common_mat_lst <- lapply(celline_group_lst, function(brcd_lst){
  ## https://stackoverflow.com/questions/1719447/outer-equivalent-for-non-vector-lists-in-r
  common_brcds <- Vectorize(function(x,y) {
    vec1 <- brcd_lst[[x]]
    vec2 <- brcd_lst[[y]]
    length(intersect(vec1,vec2))
  })
  ## Number of common barcodes across comparison-cell line in the same barcode group
  mat <- outer(1:length(brcd_lst), 1:length(brcd_lst), common_brcds)
  colnames(mat) <- names(brcd_lst)
  rownames(mat) <- names(brcd_lst)
  mat
})

pdf(glue("{plots_dir}/barcode_group_common_heatmap_ref_9.11.pdf"))
lapply(names(common_mat_lst), function(grp){
  mat <- common_mat_lst[[grp]]
  colnames(mat) <- colnames(mat) %>% str_replace("\\..+", "")
  rownames(mat) <- rownames(mat) %>% str_replace("\\..+", "")
  hm <- ComplexHeatmap::Heatmap(mat,
                                name = glue("#common {grp} barcodes"),
                                col = hcl.colors(12, "BluYl"),
                                cluster_columns=TRUE,
                                cluster_rows = TRUE,
                                row_names_gp = gpar( fontsize = 7),
                                column_names_gp = gpar(fontsize = 7),
                                show_row_names=TRUE,
                                show_column_names=TRUE,
                                #top_annotation = annot_col,
                                #left_annotation = annot_row,
                                heatmap_legend_param = list(direction = "horizontal"),
                                cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                                  grid.text(round(mat[i, j],2), x, y,gp = gpar(fontsize = 9))}
  )
  ComplexHeatmap::draw(hm,
                       merge_legend = TRUE,
                       heatmap_legend_side = "top",
                       annotation_legend_side = "top")
})
dev.off()

