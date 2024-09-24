# importing only key functions that are actually used - not to polute namespace!
import::from(.from = readr, read_csv, cols)
import::from(magrittr, "%>%")
import::from(dplyr, mutate, select, filter, rename, arrange, desc, group_by, summarise, ungroup)  # dplyr_mutate = mutate
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


experiment_name="nk_tum_immunoedit_complete" # this can be a subset analysis - e.g. just batch1,...
experiment_design="1"  # used only to construct initial dds object
abs_filt_samples=3

# parameters for annotation
biomart_host = "http://nov2020.archive.ensembl.org"
#biomart_host = "https://www.ensembl.org"
biomart_Ens_version = "Ensembl Genes 102"
biomart_dataset="mmusculus_gene_ensembl"

# directories
output_dir = paste0("results_dir", "/",experiment_name,"/")
dir.create(output_dir)





metadata <- readr::read_tsv(file = "~/workspace/datasets/metadata/quantseq_metadata_exp9.4-8.txt") 
#tx2gene_file <- file.path( "~/workspace/datasets/annotations/EnsDb.Mmusculus.v102_tx2gene.tsv")


gtf_file <- file.path("~/workspace/datasets/annotations/Mus_musculus.GRCm38.102.gtf") # needed only once to generate tx2gene.tsv
tx2gene_file <- file.path( "~/workspace/datasets/annotations/EnsDb.Mmusculus.v102_tx2gene.tsv")
tx2gene_file_md5sum <- file.path(annotation_dir, "EnsDb.Mmusculus.v102_tx2gene.md5sum")

if (!file.exists(tx2gene_file)) {
  message(tx2gene_file, " does not exist: generating!")
  txDB <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format="gtf")

  # AnnotationDbi::keys
  tx_name <- annot_db_keys(txDB, keytype="TXNAME")

  # AnnotationDbi::select
  tx2gene <- annot_db_select(txDB,
                                 keys = tx_name,
                                 columns="GENEID",
                                 keytype = "TXNAME")


  # saving tx2gene into destdir
  write.table(tx2gene,
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
  
  
  quant_files <- paste0("~/workspace/datasets/salmonquant/", metadata$processed_data_file)
  
  names(quant_files) <- metadata$library_name
  
  txi_object <- tximport::tximport(quant_files,
                                   type="salmon",
                                   tx2gene = tx2gene,
                                   ignoreTxVersion = TRUE)
  rownames(metadata) <- metadata$library_name
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

#transformation after filtering
if(!file.exists(precomputed_transf_objects_file)){
  # transformation
  log2_norm_filt <- DESeq2::normTransform(dds_filt)
  vsd_filt <- DESeq2::vst(dds_filt, blind = TRUE) # blind = TRUE for QC
  #rld_filt <- rlog(dds_filt, blind = TRUE)
  #rld_filt - too big for >100 samples; skipping!
  #rld_filt <- rlog(dds_filt, blind = FALSE) # not blind to batch effects
  save(log2_norm_filt, vsd_filt,   
       file = precomputed_transf_objects_file)
} else{
  load(precomputed_transf_objects_file)
}
