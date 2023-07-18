#' generateEnsemblAnnotation: Generates ensembl annotation file for genes in RNA-seq experiment
#' 
#' \code{generateEnsemblAnnotation} constructs, for example, DESeq2 dds object from Salmon count files, condition table and design
#' @param ensembl_ids ensembl_ids of genes in an experiment
#' @param host ensembl website address (e.g. http://www.ensembl.org)
#' @param dataset e.g. mmusculus_gene_ensembl
#' @param version version of dataset

generateEnsemblAnnotation <- function(ensembl_ids=NULL, host="http://www.ensembl.org", dataset="mmusculus_gene_ensembl", version="Ensembl Genes 92"){
  #TODO:
  # [ ] - only for mouse at the moment, add option for human (e.g. change mgi_symbol to hgnc_symbol etc.)
  # [ ] - change! by default it assumes filter type ensembl_id
  require("biomaRt") # Gene annotation using ensembl database; carafully choose same version that was used for alignment
  require("RCurl") # proxy settings for biomaRt  
  #options(RCurlOptions = list(proxy="specify-proxy-address",http.version=HTTP_VERSION_1_0)) 
  require("dplyr")
  
  ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                    dataset=dataset, 
                    host=host, 
                    port = 80, 
                    verbose =T, 
                    version=version)
  
  mart <- useDataset(dataset = dataset, ensembl)
  
  # check the available "filters" - things you can filter for
  #listFilters(ensembl) %>% filter(str_detect(name, "ensembl"))
  filterType <- "ensembl_gene_id"
  # check the available "attributes" - things you can retreive
  #listAttributes(ensembl) %>%  head(20)
  # Set the list of attributes
  
  # temporary fix to annotate also human!
  if(dataset=="hsapiens_gene_ensembl") {
    cat("Annotating: ", dataset, "\n")
    # entrezgene
    attributeNames <- c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol', 'description', 'gene_biotype', 
                        'chromosome_name', 'start_position', 'end_position', 'strand')
  } else {
    attributeNames <- c('ensembl_gene_id', 'entrezgene_id', 'mgi_symbol', 'mgi_description', 'gene_biotype',
                        'chromosome_name', 'start_position', 'end_position', 'strand')
  }

  if (!is.null(ensembl_ids) & is.character(ensembl_ids)) {
    print(paste0("Generating annotation file for genes: ", length(ensembl_ids)))
    print(paste0("dataset: ", dataset))
    print(paste0("version: ", version))
    print(paste0("host: ", host))
    
    # Adding gene names position, etc.
    ensembl_ids <- ensembl_ids
    
    # uniprot swissprot
    # add attribute: mgi_symbol
    ensemblAnnot <- getBM(
      filters = filterType, 
      attributes = attributeNames,
      values = ensembl_ids,
      mart = mart)
    
    print("Removing duplicates from ensemblAnnot file...this may take a moment!")
    ensemblAnnot_dupRemove <- ensemblAnnot %>%
      dplyr::rename(ensembl_id = ensembl_gene_id) %>%
      dplyr::group_by(ensembl_id) %>%
      dplyr::summarise_all(., .funs = function(x) {paste(unique(x), collapse = ";")}) # removing duplicates in ensembl_id by grouping; collapsing unique by ";" 
    ensemblAnnot_dupRemove <- as.data.frame(ensemblAnnot_dupRemove)
    rownames(ensemblAnnot_dupRemove) <- ensemblAnnot_dupRemove$ensembl_id
    
    # converting "NA" strings to NA
    ensemblAnnot_dupRemove$entrezgene <- ifelse(ensemblAnnot_dupRemove$entrezgene == "NA", 
                                                NA_character_, 
                                                ensemblAnnot_dupRemove$entrezgene)
    
    return(ensemblAnnot_dupRemove)
    
  } else {
    print("ensembl_gene_id needs to be specified as a character vector.")
  }
}

# other option (applicable offline and potentially faster) is to use pre-build ensembldb; 
#  Note: using ensembldb return entrez_ids as list; also it returns gene_names and not mgi_symbols (these can differ especially for older ensembl versions!)
# 
# {
#   
#   # Generate annotation file using EnsemblDB
#   test_ids <- rownames(dds)#[1:10]                                                                                                                                                                      
#   
#   ## We're going to fetch all genes which names start with BCL.   
#   require(EnsDb.Mmusculus.v95)
#   test_ids_annot <- genes(EnsDb.Mmusculus.v95,                                                                                                                                                                          
#                           columns = c("gene_id", "gene_name","entrezid", "gene_biotype", "description"), #"gene_name", "symbol",  "entrezid",                                                                                               
#                           filter = GeneIdFilter(test_ids, condition = "=="),                                                                                                                                                     
#                           return.type = "data.frame") # return.type = "DataFrame"; DataFrame-class {S4Vectors}                                                                                                                   
#   test_ids_annot_df <- as.data.frame(test_ids_annot)                                                                                                                                                                   
#   
#   require(biomaRt)
#   filters_param = 'ensembl_gene_id'
#   attrib_param = c('ensembl_gene_id', 'mgi_symbol', 'entrezgene','gene_biotype','description') # description - associated name description ENSEMBL                                                                                   
#   ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset=biomart_dataset,
#                     host=biomart_host,port = 80, verbose = T, 
#                     version = biomart_Ens_version)  
#   
#   tgondii_txdb <- GenomicFeatures::makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
#                                                        dataset=biomart_dataset,
#                                                        host=biomart_host,port = 80)
#   
#   annotation_df <- getBM(attributes=attrib_param, filters=filters_param, values=test_ids, mart=ensembl) #, verbose = T
#   annotation_df_duplicated <- annotation_df[duplicated(annotation_df$ensembl_gene_id) | duplicated(annotation_df$ensembl_gene_id, fromLast=TRUE),]
#   test_ids_annot2 <- test_ids_annot                                                                                                                                                                                    
#   colnames(test_ids_annot2) <- c("ensembl_gene_id", "mgi_symbol", "gene_biotype", "description")   
#   
#   # genes() return gene_names not mgi_symbols; but latest ensembl seems to mainly use mgi_symbols; but for an archive version gene_name seems to be correct
#   # getBM returns 
#   
#   dplyr::all_equal(test_ids_annot2, annotation_df)
#   identical(test_ids_annot2, annotation_df) # not identical, because there are duplicates in annotation_df                                                                                                             
#   all.equal(test_ids_annot2, annotation_df)                    
# }
